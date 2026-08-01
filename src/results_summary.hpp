// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_RESULTS_SUMMARY_HPP_
#define SRC_RESULTS_SUMMARY_HPP_

#include "adapter_matcher.hpp"
#include "base_groups.hpp"
#include "duplication_results.hpp"
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "file_info.hpp"
#include "kmer_counter.hpp"
#include "quality_score.hpp"
#include "results_collector.hpp"
#include "run_mode.hpp"
#include "tile_processor.hpp"

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

struct results_summary {
  std::uint64_t n_reads{};
  std::uint64_t min_read_len{};
  std::uint64_t max_read_len{};
  std::uint64_t median_read_len{};
  std::uint64_t total_bases{};
  std::uint64_t total_gc{};
  std::vector<falco::nuc_array> base_counts;
  falco::gc_content_array gc_content{};
  std::vector<std::uint64_t> n_counts;
  std::vector<std::uint64_t> lengths;
  std::vector<falco::qual_array> qual_by_pos;
  falco::qual_array qual_by_read{};
  duplication_results dr;
  std::uint64_t n_reads_for_dups{};
  std::vector<overrep_t> overrep;
  dup_summary_t dup_summary;
  adapter_matcher am;
  tile_processor tp;
  tile_processor::tiles_centered_t centered;
  kmer_counter kc;
  std::vector<kmer_result> kmer_results;
  file_grades grades;
  run_mode mode;
  file_info info;
  base_group_vec groups;

  // clang-format off
  results_summary(results_collector &rc, const run_mode &mode,
                  const file_info &info) :
    n_reads{rc.n_reads},
    max_read_len{rc.max_read_len},
    base_counts{std::move(rc.base_counts)},
    gc_content{rc.gc_content},
    n_counts{std::move(rc.n_counts)},
    lengths{std::move(rc.lengths)},
    qual_by_pos{std::move(rc.qual_by_pos)},
    qual_by_read{rc.qual_by_read},
    dr{std::move(rc.dr)},
    am{std::move(rc.am)},
    tp{std::move(rc.tp)},
    kc{std::move(rc.kc)},
    mode{mode},
    info{info} {
    initialize();
  }

  auto
  apply_groups() -> void;

  auto
  assign_grades() -> void;

  auto
  initialize() -> void;

  [[nodiscard]] auto
  get_report() const -> std::string;

  [[nodiscard]] auto
  get_html() const -> std::string;

  [[nodiscard]] auto
  get_summary() const -> std::string;
};

#endif  // SRC_RESULTS_SUMMARY_HPP_
