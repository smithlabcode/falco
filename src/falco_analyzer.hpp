// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_ANALYZER_HPP_
#define SRC_FALCO_ANALYZER_HPP_

#include "results_collector.hpp"
#include "task_queue.hpp"

#include <atomic>
#include <cstdint>
#include <format>
#include <vector>

class reads_file_t;
class run_mode;
struct file_info;

struct analyzer_t {
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
