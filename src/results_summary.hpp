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

#ifndef SRC_RESULTS_SUMMARY_HPP_
#define SRC_RESULTS_SUMMARY_HPP_

#include "adapter_matcher.hpp"
#include "duplication_results.hpp"
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "kmer_counter.hpp"
#include "quality_score.hpp"
#include "results_collector.hpp"
#include "run_mode.hpp"
#include "tile_processor.hpp"

#include <cstdint>
#include <format>
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
  std::vector<base_group_t> groups;

  results_summary(results_collector &rc, const run_mode &mode,
                  const file_info &info) : mode{mode}, info{info} {
    n_reads = rc.n_reads;
    max_read_len = rc.max_read_len;
    base_counts = std::move(rc.base_counts);
    gc_content = std::move(rc.gc_content);
    n_counts = std::move(rc.n_counts);
    lengths = std::move(rc.lengths);
    qual_by_pos = std::move(rc.qual_by_pos);
    qual_by_read = std::move(rc.qual_by_read);
    dr = std::move(rc.dr);
    am = std::move(rc.am);
    initialize();
  }

  results_summary(results_collector_tile &rc, const run_mode &mode,
                  const file_info &info) : mode{mode}, info{info} {
    n_reads = rc.n_reads;
    max_read_len = rc.max_read_len;
    base_counts = std::move(rc.base_counts);
    gc_content = std::move(rc.gc_content);
    n_counts = std::move(rc.n_counts);
    lengths = std::move(rc.lengths);
    qual_by_pos = std::move(rc.qual_by_pos);
    qual_by_read = std::move(rc.qual_by_read);
    dr = std::move(rc.dr);
    am = std::move(rc.am);
    //
    tp = std::move(rc.tp);

    initialize();
  }

  results_summary(results_collector_kmer &rc, const run_mode &mode,
                  const file_info &info) : mode{mode}, info{info} {
    n_reads = rc.n_reads;
    max_read_len = rc.max_read_len;
    base_counts = std::move(rc.base_counts);
    gc_content = std::move(rc.gc_content);
    n_counts = std::move(rc.n_counts);
    lengths = std::move(rc.lengths);
    qual_by_pos = std::move(rc.qual_by_pos);
    qual_by_read = std::move(rc.qual_by_read);
    dr = std::move(rc.dr);
    am = std::move(rc.am);
    //
    kc = std::move(rc.kc);

    initialize();
  }

  results_summary(results_collector_tile_kmer &rc, const run_mode &mode,
                  const file_info &info) : mode{mode}, info{info} {
    n_reads = rc.n_reads;
    max_read_len = rc.max_read_len;
    base_counts = std::move(rc.base_counts);
    gc_content = std::move(rc.gc_content);
    n_counts = std::move(rc.n_counts);
    lengths = std::move(rc.lengths);
    qual_by_pos = std::move(rc.qual_by_pos);
    qual_by_read = std::move(rc.qual_by_read);
    dr = std::move(rc.dr);
    am = std::move(rc.am);
    //
    tp = std::move(rc.tp);
    kc = std::move(rc.kc);

    initialize();
  }

  auto
  apply_groups(const std::vector<base_group_t> &groups) -> void;

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
