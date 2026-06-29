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

#ifndef SRC_DUPLICATION_RESULTS_HPP_
#define SRC_DUPLICATION_RESULTS_HPP_

#include "falco_word.hpp"

#include "boost/boost_unordered.hpp"

#include <cstdint>
#include <string>
#include <vector>

class run_mode;
struct file_grades;
struct file_info;

struct overrep_t {
  falco_word w;
  std::uint64_t n_obs{};     // number of observations
  double pct_val{};          // corresponding percentage
  std::int32_t contam_id{};  // index of matching contaminant
};

struct dup_summary_t {
  std::uint64_t max_dup{};
  std::uint64_t n_reads{};
  std::vector<std::uint64_t> hist_mass;
  std::vector<std::uint64_t> hist_dedup;
};

struct duplication_results {
  static constexpr auto max_n_reads_total{1'000'000};
#ifdef ORIGINAL_DUPS
  static constexpr auto max_reads_to_hash{100'000};
#endif  // ORIGINAL_DUPS
  static constexpr auto default_read_skip{10};
  static constexpr auto overrep_cutoff = 0.001;

  std::int64_t read_skip{default_read_skip};
  std::int64_t read_idx{};
  boost::unordered_flat_map<falco_word, std::uint64_t> dups;

  auto
  initialize(const run_mode &mode, const file_info &info) -> void;

  auto
  get_n_reads() const -> std::uint64_t;

  [[nodiscard]] auto
  get_dups_summary() const -> dup_summary_t;

  [[nodiscard]] auto
  get_overrepresented(const std::uint64_t n_reads) const
    -> std::vector<overrep_t>;

  auto
  operator+=(const duplication_results &rhs) -> const duplication_results &;

#ifdef ORIGINAL_DUPS
  auto
  count_seqs(const auto seq_itr, const auto sz) {
    const auto fw = falco_word(seq_itr, sz);
    if (read_idx++ < max_reads_to_hash || dups.contains(fw))
      ++dups[fw];
  }
#else   // NOT ORIGINAL_DUPS
  auto
  count_seqs(const auto seq_itr, const auto sz) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      ++dups[falco_word(seq_itr, sz)];
    }
  }
#endif  // ORIGINAL_DUPS
};

[[nodiscard]] auto
get_grade_duplication(const dup_summary_t &summary) -> std::string;

[[nodiscard]] auto
duplication_report(const dup_summary_t &summary,
                   const file_grades &grades) -> std::string;

[[nodiscard]] auto
duplication_html(const dup_summary_t &summary,
                 const file_grades &grades) -> std::string;

[[nodiscard]] auto
get_grade_overrepresented(const std::uint64_t n_reads,
                          const duplication_results &dr) -> std::string;

[[nodiscard]] auto
overrepresented_report(const std::vector<overrep_t> &overrep,
                       const file_grades &grades) -> std::string;

[[nodiscard]] auto
overrepresented_html(const std::vector<overrep_t> &overrep,
                     const file_grades &grades) -> std::string;

#endif  // SRC_DUPLICATION_RESULTS_HPP_
