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
  std::mutex queue_mtx;
  std::condition_variable_any queue_ready;  // can be canceled
  std::vector<std::mutex> file_mtx;
  std::condition_variable file_ready;
  std::uint32_t queue_size{};
  std::queue<task_t> tasks;
  std::vector<std::int32_t> n_tasks;
  std::vector<std::vector<results_t>> results;

  explicit analyzer_t(const std::uint32_t n_threads,
                      const std::uint32_t n_files, auto &reads_files,
                      const std::vector<file_info> &infos) :
    file_mtx(n_files), queue_size(n_chunks_per_thread * n_threads * n_files),
    n_tasks(n_files), results(n_threads, std::vector<results_t>(n_files)) {
    assert(std::size(reads_files) == std::size(infos));

    // set any per-file information used to do the analysis
    for (auto &res : results)
      for (const auto [file_id, info] : std::views::enumerate(infos))
        res[file_id].init(info);

    const auto n_chunks_per_file = n_threads * n_chunks_per_thread;

    std::vector<std::jthread> workers;
    workers.reserve(n_threads);
    for (auto th_id = 0u; th_id < n_threads; ++th_id) {
      workers.emplace_back([th_id, this](const std::stop_token &stop) {
        auto &r = results[th_id];
        task_t task;
        while (true) {
          std::unique_lock ql(queue_mtx);
          queue_ready.wait(ql, stop, [this] { return !tasks.empty(); });
          if (stop.stop_requested() && tasks.empty())
            return;
          task = std::move(tasks.front());
          tasks.pop();
          ql.unlock();

          // both readers and workers might be waiting on this
          queue_ready.notify_all();

          r[task.file_id].template process_reads<rec_t>(task.beg, task.end);

          std::unique_lock fl(file_mtx[task.file_id]);
          --n_tasks[task.file_id];
          fl.unlock();

          // for readers waiting to delete owned data; notify all because we
          // don't know which one
          file_ready.notify_all();
        }
      });
    }

    std::vector<std::thread> readers;
    readers.reserve(n_files);
    for (auto f_id = 0u; f_id < n_files; ++f_id)
      readers.emplace_back([&, f_id]() {
        while (reads_files[f_id]) {
          std::unique_lock fl(file_mtx[f_id]);
          file_ready.wait(fl, [f_id, this] { return n_tasks[f_id] == 0; });
          fl.unlock();
          reads_files[f_id].load_next();
          const auto chunks = get_chunks(reads_files[f_id], n_chunks_per_file);
          std::ranges::for_each(chunks,
                                [&](const auto &c) { push_task(f_id, c); });
        }
      });
    std::ranges::for_each(readers, [](auto &r) { r.join(); });

    // tell workers to stop
    std::unique_lock ql(queue_mtx);
    queue_ready.wait(ql, [this] { return tasks.empty(); });
    std::ranges::for_each(workers, [](auto &w) { w.request_stop(); });
    ql.unlock();
    queue_ready.notify_all();
  }

  auto
  push_task(const auto file_id, const auto &chunk) {
    // wait until the queue has space
    std::unique_lock ql(queue_mtx);
    queue_ready.wait(ql, [this] { return std::size(tasks) < queue_size; });
    tasks.emplace(file_id, chunk.first, chunk.second);
    ql.unlock();
    queue_ready.notify_one();
    std::lock_guard fl(file_mtx[file_id]);
    ++n_tasks[file_id];  // no notify: only this thread waits on this value

    // any number of workers could be waiting to update n_tasks[file_id]
    file_ready.notify_all();
  }
};

#endif  // SRC_FALCO_ANALYZER_HPP_
