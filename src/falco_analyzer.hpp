// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_ANALYZER_HPP_
#define SRC_FALCO_ANALYZER_HPP_

#include "file_info.hpp"
#include "reads_file.hpp"
#include "results_collector.hpp"
#include "run_mode.hpp"
#include "task_queue.hpp"

#include <atomic>
#include <cstdint>
#include <vector>

struct analyzer_t {
  static constexpr auto n_chunks_per_thread = 8;  // ADS: not empirically tested
  task_queue tq;
  // ADS: likely false sharing between n_tasks members
  std::vector<std::atomic_int32_t> n_tasks;  // count of submitted tasks by file
  std::atomic_uint32_t n_active_files{};
  std::vector<std::vector<results_collector>> results;  // n_threads x n_tiles
  analyzer_t(const std::uint32_t n_threads,             //
             const run_mode &mode,                      //
             const std::vector<file_info> &infos,       //
             std::vector<reads_file_t> &reads_files     //
  );
};

#endif  // SRC_FALCO_ANALYZER_HPP_
