// SPDX-License-Identifier: MIT; (c) 2026 Andrew D Smith

#ifndef SRC_FALCO_ANALYZER_HPP_
#define SRC_FALCO_ANALYZER_HPP_

#include "bam_file.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "file_info.hpp"
#include "reads_file.hpp"
#include "run_mode.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <condition_variable>
#include <cstdint>
#include <format>
#include <iterator>
#include <list>
#include <map>
#include <mutex>
#include <ranges>
#include <stop_token>
#include <string>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

template <typename results_t> struct analyzer_t {
  static constexpr auto n_chunks_per_thread = 8;  // ADS: not empirically tested

  std::mutex mtx;
  std::condition_variable cv;
  std::list<std::pair<std::int32_t, task_t>> q;  // file_id and task

  std::vector<std::atomic_int32_t> n_tasks;  // count of submitted tasks by file
  std::uint32_t n_active_files{};
  std::vector<std::vector<results_t>> results;
  bool stop{};

  analyzer_t(const std::uint32_t n_threads,          //
             const std::uint32_t n_files,            //
             const run_mode &mode,                   //
             const std::vector<file_info> &infos,    //
             std::vector<reads_file_t> &reads_files  //
  );

  auto
  push(const auto file_id, task_t task) {
    std::unique_lock tq_lock(mtx);
    q.emplace_back(file_id, std::move(task));
    tq_lock.unlock();
    cv.notify_one();
  }

  auto
  push_each(const auto f_id, auto tasks) {
    std::ranges::for_each(tasks, [&](auto &&c) { push(f_id, std::move(c)); });
  }

  [[nodiscard]] auto
  pop_and_release_lock(auto &tq_lock) {
    auto task = std::move(q.front());
    q.pop_front();
    tq_lock.unlock();
    cv.notify_one();
    return task;
  }
};

template <typename results_t>
analyzer_t<results_t>::analyzer_t(const std::uint32_t n_threads,
                                  const std::uint32_t n_files,
                                  const run_mode &mode,
                                  const std::vector<file_info> &infos,
                                  std::vector<reads_file_t> &reads_files) :
  n_tasks(n_files), n_active_files{n_files},
  results(n_threads, std::vector<results_t>(n_files)) {
  assert(std::size(reads_files) == std::size(infos));
  const auto n_chunks_per_file = n_threads * n_chunks_per_thread;
  // set per-file information used to do the analysis
  for (auto &res : results)
    for (const auto [file_id, info] : std::views::enumerate(infos))
      res[file_id].init(mode, info);

  // add initial jobs so workers can work immediately
  for (const auto f_id : std::views::iota(0u, n_files))
    push(f_id, std::monostate{});

  std::vector<std::jthread> workers;
  workers.reserve(n_threads);
  for (const auto th_id : std::views::iota(0u, n_threads))
    workers.emplace_back([this, n_chunks_per_file, th_id, &reads_files] {
      auto &res = results[th_id];
      while (true) {
        std::unique_lock tq_lock(mtx);
        cv.wait(tq_lock, [this] { return !q.empty() || stop; });
        if (stop && q.empty())
          return;
        auto [file_id, task] = pop_and_release_lock(tq_lock);
        if (std::holds_alternative<fq_task_t>(task))
          process_reads_fq(res[file_id], std::get<fq_task_t>(task));
        else if (std::holds_alternative<bam_task_t>(task))
          process_reads_bam(res[file_id], std::get<bam_task_t>(task));
        else if (std::holds_alternative<bgzf_block_t>(task))
          decompress(std::get<bgzf_block_t>(task));
        else {  // monostate means read more data
          if (!is_active(reads_files[file_id])) {
            if (--n_active_files == 0) {
              stop = true;
              cv.notify_all();
              return;
            }
          }
          else {
            if (inflate_only(reads_files[file_id])) {
              auto tasks = make_tasks_inflate(reads_files[file_id]);
              n_tasks[file_id] = std::ssize(tasks) + 1;  // +1 for self
              push_each(file_id, std::move(tasks));
            }
            else {
              auto tasks = make_tasks(reads_files[file_id], n_chunks_per_file);
              n_tasks[file_id] = std::ssize(tasks) + 1;  // +1 for self
              push_each(file_id, std::move(tasks));
              tasks = make_tasks_inflate(reads_files[file_id]);
              n_tasks[file_id] += std::ssize(tasks);
              push_each(file_id, std::move(tasks));
            }
          }
        }
        if (n_tasks[file_id].fetch_sub(1) == 1)
          push(file_id, std::monostate{});
      }
    });
}

#endif  // SRC_FALCO_ANALYZER_HPP_
