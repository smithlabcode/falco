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
#include "contaminants.hpp"
#include "falco_file_format.hpp"
#include "falco_results.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "quality_score.hpp"
#include "tile_processor.hpp"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <print>
#include <queue>
#include <ranges>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

template <typename results_t, typename rec_t> struct analyzer_t {
  static constexpr auto n_chunks_per_thread = 8;
  using pos_t = typename rec_t::pos_t;
  struct task_t {
    std::int64_t file_id{};
    pos_t beg{};
    pos_t end{};
  };
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
                      const std::uint32_t n_files, auto &reads_files,
                      const std::vector<file_info> &infos) :
    n_tasks(n_files), n_active_files{n_files},
    results(n_threads, std::vector<results_t>(n_files)) {
    assert(std::size(reads_files) == std::size(infos));
    const auto n_chunks_per_file = n_threads * n_chunks_per_thread;
    // set per-file information used to do the analysis
    for (auto &res : results)
      for (const auto [file_id, info] : std::views::enumerate(infos))
        res[file_id].init(info);

    for (const auto f_id : std::views::iota(0u, n_files))
      file_queue.emplace(f_id);

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
          const auto task = pop_task(task_queue_lock);
          res[task.file_id].template process_reads<rec_t>(task.beg, task.end);
          if (n_tasks[task.file_id].fetch_sub(1) == 1)
            push_file(task.file_id);
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
          const auto f_id = pop_file(file_queue_lock);
          if (!reads_files[f_id]) {
            if (--n_active_files == 0) {
              file_queue_cv.notify_all();
              return;
            }
            continue;
          }
          reads_files[f_id].load_next();
          n_tasks[f_id] = n_chunks_per_file;
          const auto chunks = get_chunks(reads_files[f_id], n_chunks_per_file);
          std::ranges::for_each(chunks,
                                [&](const auto &c) { push_task(f_id, c); });
        }
      });
    std::ranges::for_each(readers, [](auto &reader) { reader.join(); });

    std::unique_lock task_queue_lock(task_queue_mtx);
    task_queue_cv.wait(task_queue_lock, [this] { return task_queue.empty(); });
    std::ranges::for_each(workers, [](auto &worker) { worker.request_stop(); });
  }

  auto
  push_task(const auto file_id, const auto &chunk) {
    std::unique_lock task_queue_lock(task_queue_mtx);
    task_queue.emplace(file_id, chunk.first, chunk.second);
    task_queue_lock.unlock();
    task_queue_cv.notify_one();
  }

  [[nodiscard]] auto
  pop_task(auto &task_queue_lock) {
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
  pop_file(auto &file_queue_lock) {
    const auto file_id = std::move(file_queue.front());
    file_queue.pop();
    file_queue_lock.unlock();
    file_queue_cv.notify_one();
    return file_id;
  }
};

#endif  // SRC_FALCO_ANALYZER_HPP_
