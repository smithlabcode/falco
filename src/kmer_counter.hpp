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

#ifndef SRC_KMER_COUNTER_HPP_
#define SRC_KMER_COUNTER_HPP_

#include "falco_utils.hpp"

#include <array>
#include <cstdint>
#include <format>
#include <iterator>
#include <limits>
#include <string>
#include <thread>  // IWYU pragma: keep
#include <utility>
#include <vector>

struct kmer_result {
  // holds info for {kmer} x {position} to be sorted, filtered and output
  std::uint64_t kmer{};
  std::uint64_t count{};
  double pval{1.0};
  double obs_exp{};
  std::uint64_t pos{};
  [[nodiscard]] auto
  operator<=>(const kmer_result &rhs) const;
  [[nodiscard]] auto
  string() const;
};

struct kmer_counter {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{2.0, "pass"},
    std::pair{5.0, "warn"},
    std::pair{std::numeric_limits<double>::max(), "fail"},
  };

  // cutoffs for what to report
  static constexpr auto min_obs_exp_to_report = 5.0;
  static constexpr auto max_pval_to_report = 0.01;
  [[maybe_unused]] static constexpr auto max_kmers_to_plot = 10;
  static constexpr auto n_kmers_to_report = 20;

  // don't count kmers past this length of read prefix
  static constexpr auto max_len_for_analysis = 500;

  static constexpr auto read_skip = 50 - 1;
  static constexpr auto kmer_size = 7;
  static constexpr auto n_kmers = ipow(4, kmer_size);
  static constexpr auto kmer_mask = n_kmers - 1;

  std::int32_t read_idx{};
  std::int32_t max_read_len{};
  std::vector<std::array<std::uint64_t, n_kmers>> kmer_counts;

  auto
  resize(const std::int32_t updated_length) {
    max_read_len = updated_length;
    kmer_counts.resize(updated_length);
  }

  static constexpr auto
  shift_kmer(auto &k, auto c) {
    k = (k << 2) | encode(c);
  };

  auto
  count_kmers(auto seq_itr, auto sz) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      sz = sz > max_len_for_analysis ? max_len_for_analysis : sz;
      if (sz > max_read_len) [[unlikely]]
        resize(sz);
      const auto seq_end = seq_itr + sz;
      auto out_itr = std::begin(kmer_counts);
      auto kmer = 0;
      for (auto i = 0; i < kmer_size - 1 && seq_itr != seq_end; ++i)
        shift_kmer(kmer, *seq_itr++);
      while (seq_itr != seq_end) {
        shift_kmer(kmer, *seq_itr++);
        // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-constant-array-index)
        ++(*out_itr++)[kmer & kmer_mask];
      }
    }
  }

  auto
  operator+=(const kmer_counter &rhs) -> const kmer_counter &;

  [[nodiscard]] auto
  get_kmer_results() const -> std::vector<kmer_result>;

  [[nodiscard]] auto
  get_report(std::string &grade) const -> std::string;

  [[nodiscard]] static auto
  decode_kmer(auto word, const auto n_bases) -> std::string;
};

#endif  // SRC_KMER_COUNTER_HPP_
