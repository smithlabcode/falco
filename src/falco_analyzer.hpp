// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_ANALYZER_HPP_
#define SRC_FALCO_ANALYZER_HPP_

#include "bam_file.hpp"
#include "falco_task.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "file_info.hpp"
#include "reads_file.hpp"
#include "run_mode.hpp"
#include "task_queue.hpp"

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

  task_queue tq;

  std::vector<std::atomic_int32_t> n_tasks;  // count of submitted tasks by file
  std::uint32_t n_active_files{};
  std::vector<std::vector<results_t>> results;

  analyzer_t(const std::uint32_t n_threads,          //
             const std::uint32_t n_files,            //
             const run_mode &mode,                   //
             const std::vector<file_info> &infos,    //
             std::vector<reads_file_t> &reads_files  //
  );
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
  const auto n_chunks = n_threads * n_chunks_per_thread;
  // set per-file information used to do the analysis
  for (auto &res : results)
    for (const auto [file_id, info] : std::views::enumerate(infos))
      res[file_id].init(mode, info);

  // add initial jobs so workers can work immediately
  for (const auto file_id : std::views::iota(0u, n_files))
    tq.push(file_id, std::monostate{});

  std::vector<std::jthread> workers;
  workers.reserve(n_threads);
  for (const auto th_id : std::views::iota(0u, n_threads))
    workers.emplace_back([this, n_chunks, th_id, &reads_files] {
      auto &res = results[th_id];
      while (true) {
        auto tq_lock = tq.wait_and_acquire_lock();
        if (tq.is_finished())
          return;
        auto [file_id, task] = tq.pop_and_release_lock(tq_lock);
        if (std::holds_alternative<fq_task_t>(task))
          process_reads_fq(res[file_id], std::get<fq_task_t>(task));
        else if (std::holds_alternative<bam_task_t>(task))
          process_reads_bam(res[file_id], std::get<bam_task_t>(task));
        else if (std::holds_alternative<bgzf_block_t>(task))
          decompress(std::get<bgzf_block_t>(task));
        else {  // monostate means read more data
          if (!is_active(reads_files[file_id])) {
            if (--n_active_files == 0) {
              tq.request_shutdown();
              return;
            }
          }
          else
            make_tasks(reads_files[file_id], n_chunks, file_id, tq,
                       n_tasks[file_id]);
        }
        if (n_tasks[file_id].fetch_sub(1) == 1)
          tq.push(file_id, std::monostate{});
      }
    });
}

#endif  // SRC_FALCO_ANALYZER_HPP_
