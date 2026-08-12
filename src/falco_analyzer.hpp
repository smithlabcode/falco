// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_ANALYZER_HPP_
#define SRC_FALCO_ANALYZER_HPP_

#include "results_collector.hpp"

#include <cstdint>
#include <vector>

class reads_file_t;
class run_mode;
struct dups_init_t;
struct file_info;

[[nodiscard]] auto
analyze(const std::uint32_t n_threads, const run_mode &mode,
        std::vector<file_info> &infos, std::vector<reads_file_t> reads_files,
        std::vector<dups_init_t> dups_init) -> std::vector<results_collector>;

#endif  // SRC_FALCO_ANALYZER_HPP_
