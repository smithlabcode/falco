// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_GC_SUMMARY_HPP_
#define SRC_FALCO_GC_SUMMARY_HPP_

#include "falco_utils.hpp"

#include <cstdint>
#include <vector>

[[nodiscard]] auto
combine_gc_content_for_lengths(const std::vector<falco::gc_content_t> &gcs)
  -> std::vector<double>;

[[nodiscard]] auto
smooth_gc_content(const std::vector<double> &data,
                  const std::int64_t window_size) -> std::vector<double>;

[[nodiscard]] auto
get_theoretical_distribution(const std::vector<double> &gc,
                             const std::uint64_t total_count)
  -> std::vector<double>;

[[nodiscard]] auto
sum_deviation_from_normal(const std::vector<double> &gc) -> double;

#endif  // SRC_FALCO_GC_SUMMARY_HPP_
