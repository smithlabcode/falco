/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SRC_FALCO_ANALYZER_HPP_
#define SRC_FALCO_ANALYZER_HPP_

#include "bam_file.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "file_info.hpp"
#include "run_mode.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <cstdint>
#include <format>
#include <iterator>
#include <map>
#include <mutex>
#include <queue>
#include <ranges>
#include <stop_token>
#include <string>
#include <thread>
#include <utility>
#include <vector>

struct fq_task_t {
  std::int64_t file_id{};
  fqrec::pos_t beg{};
  fqrec::pos_t end{};
};
struct bam_task_t {
  std::int64_t file_id{};
  bamrec::pos_t beg{};
  bamrec::pos_t end{};
};

using task_t = std::variant<fq_task_t, bam_task_t>;

static auto
make_tasks_impl(const std::int64_t file_id,
                const bam_chunks_t &chunks) -> std::vector<task_t> {
  std::vector<task_t> tasks;
  tasks.reserve(std::size(chunks));
  for (const auto &chunk : chunks)
    tasks.push_back(bam_task_t(file_id, chunk.first, chunk.second));
  return tasks;
}

[[nodiscard]] static auto
make_tasks_impl(const std::int64_t file_id,
                const fq_chunks_t &chunks) -> std::vector<task_t> {
  std::vector<task_t> tasks;
  tasks.reserve(std::size(chunks));
  for (const auto &chunk : chunks)
    tasks.push_back(fq_task_t(file_id, chunk.first, chunk.second));
  return tasks;
}

template <typename T>
concept concrete_reads_file_t =        //
  std::same_as<T, fastq_file> ||       //
  std::same_as<T, fastq_bgzf_file> ||  //
  std::same_as<T, fastq_gz_file> ||    //
  std::same_as<T, bam_file>;

[[nodiscard]] inline auto
make_tasks(concrete_reads_file_t auto &reads_file, const std::int64_t file_id,
           const std::int64_t n_chunks) -> std::vector<task_t> {
  reads_file.load_next();
  const auto chunks = reads_file.get_chunks(n_chunks);
  return make_tasks_impl(file_id, chunks);
}

inline auto
load_next(concrete_reads_file_t auto &f) -> void {
  f.load_next();
}

[[nodiscard]] inline auto
is_active(const concrete_reads_file_t auto &reads_file) -> bool {
  return static_cast<bool>(reads_file);
}

class reads_file_t {
public:
  template <typename T> reads_file_t(T x) : self_(new model<T>(std::move(x))) {}
  reads_file_t(const reads_file_t &x) = delete;
  auto
  operator=(const reads_file_t &x) -> reads_file_t & = delete;
  auto
  operator=(reads_file_t &&) noexcept -> reads_file_t & = delete;
  reads_file_t(reads_file_t &&) noexcept = default;
  [[nodiscard]] friend auto
  is_active(const reads_file_t &x) -> bool {
    return x.self_->is_active_();
  }
  [[nodiscard]] friend auto
  make_tasks(auto &reads_file, const std::int64_t file_id,
             const std::int64_t n_chunks) -> std::vector<task_t> {
    return reads_file.self_->make_tasks_(file_id, n_chunks);
  }

private:
  struct concept_t {
    virtual ~concept_t() = default;
    virtual auto
    make_tasks_(const std::int64_t,
                const std::int64_t) -> std::vector<task_t> = 0;
    virtual auto
    is_active_() const -> bool = 0;
  };
  template <typename T> struct model : concept_t {
    model(T x) : data_(std::move(x)) {}
    [[nodiscard]] auto
    make_tasks_(const std::int64_t file_id,
                const std::int64_t n_chunks) -> std::vector<task_t> override {
      return make_tasks(data_, file_id, n_chunks);
    }
    [[nodiscard]] auto
    is_active_() const -> bool override {
      return is_active(data_);
    }
    T data_;
  };
  std::unique_ptr<concept_t> self_;
};

template <typename results_t> struct analyzer_t {
  static constexpr auto n_chunks_per_thread = 8;
  std::mutex task_queue_mtx;
  std::condition_variable_any task_queue_cv;
  std::queue<task_t> task_queue;

  std::mutex file_queue_mtx;
  std::condition_variable file_queue_cv;
  std::queue<std::int32_t> file_queue;

  std::vector<std::atomic_int32_t> n_tasks;
  std::uint32_t n_active_files{};
  std::vector<std::vector<results_t>> results;

  explicit analyzer_t(const std::uint32_t n_threads,
                      const std::uint32_t n_readers,
                      const std::uint32_t n_files,  //
                      const run_mode &mode,
                      const std::vector<file_info> &infos,  //
                      auto &reads_files) :
    n_tasks(n_files), n_active_files{n_files},
    results(n_threads, std::vector<results_t>(n_files)) {
    assert(std::size(reads_files) == std::size(infos));
    const auto n_chunks_per_file = n_threads * n_chunks_per_thread;
    // set per-file information used to do the analysis
    for (auto &res : results)
      for (const auto [file_id, info] : std::views::enumerate(infos))
        res[file_id].init(mode, info);

    // put each file in the file queue to start
    for (const auto f_id : std::views::iota(0u, n_files))
      file_queue.emplace(f_id);

    {
      std::vector<std::jthread> workers;
      workers.reserve(n_threads);
      for (const auto th_id : std::views::iota(0u, n_threads))
        workers.emplace_back([this, th_id](const std::stop_token &stop) {
          auto &res = results[th_id];
          while (true) {
            std::unique_lock task_queue_lock(task_queue_mtx);
            task_queue_cv.wait(task_queue_lock, stop,
                               [this] { return !task_queue.empty(); });
            if (stop.stop_requested() && task_queue.empty())
              return;
            const auto task = pop_task_and_release_lock(task_queue_lock);
            auto file_id = 0L;
            if (std::holds_alternative<fq_task_t>(task)) {
              const auto fq_task = std::get<fq_task_t>(task);
              file_id = fq_task.file_id;
              res[file_id].process_reads_fq(fq_task.beg, fq_task.end);
            }
            else {  // std::holds_alternative<bam_task_t>(task)
              const auto bam_task = std::get<bam_task_t>(task);
              file_id = bam_task.file_id;
              res[file_id].process_reads_bam(bam_task.beg, bam_task.end);
            }
            if (n_tasks[file_id].fetch_sub(1) == 1)
              push_file(file_id);
          }
        });

      std::vector<std::jthread> readers;
      readers.reserve(n_readers);
      for (const auto _ : std::views::iota(0u, n_readers))
        readers.emplace_back([&]() {
          while (true) {
            std::unique_lock file_queue_lock(file_queue_mtx);
            file_queue_cv.wait(file_queue_lock, [this] {
              return n_active_files == 0 || !file_queue.empty();
            });
            if (n_active_files == 0)
              return;
            const auto file_id = pop_file_and_release_lock(file_queue_lock);
            if (!is_active(reads_files[file_id])) {
              if (--n_active_files == 0) {
                file_queue_cv.notify_all();
                return;
              }
              continue;
            }
            const auto tasks =
              make_tasks(reads_files[file_id], file_id, n_chunks_per_file);
            n_tasks[file_id] = static_cast<std::int32_t>(std::ssize(tasks));
            std::ranges::for_each(tasks, [&](const auto &c) { push_task(c); });
          }
        });
      std::ranges::for_each(readers, [](auto &reader) { reader.join(); });

      std::unique_lock task_queue_lock(task_queue_mtx);
      task_queue_cv.wait(task_queue_lock,
                         [this] { return task_queue.empty(); });
      std::ranges::for_each(workers,
                            [](auto &worker) { worker.request_stop(); });
    }
  }

  auto
  push_task(const auto &task) {
    std::unique_lock task_queue_lock(task_queue_mtx);
    task_queue.emplace(task);
    task_queue_lock.unlock();
    task_queue_cv.notify_one();
  }

  [[nodiscard]] auto
  pop_task_and_release_lock(auto &task_queue_lock) {
    const auto task = std::move(task_queue.front());
    task_queue.pop();
    task_queue_lock.unlock();
    task_queue_cv.notify_one();
    return task;
  }

  auto
  push_file(const auto file_id) {
    std::unique_lock file_queue_lock(file_queue_mtx);
    file_queue.emplace(file_id);
    file_queue_lock.unlock();
    file_queue_cv.notify_one();
  }

  [[nodiscard]] auto
  pop_file_and_release_lock(auto &file_queue_lock) {
    const auto file_id = file_queue.front();
    file_queue.pop();
    file_queue_lock.unlock();
    file_queue_cv.notify_one();
    return file_id;
  }
};

#endif  // SRC_FALCO_ANALYZER_HPP_
