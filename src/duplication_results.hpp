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

#include <array>
#include <cstdint>
#include <format>
#include <iterator>
#include <string>
#include <unordered_map>
#include <utility>

struct duplication_results {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.50, "error"},
    std::pair{0.70, "warn"},
    std::pair{1.00, "pass"},
  };

  static constexpr auto overrep_cutoffs = std::array{
    std::pair{0.01, "pass"},
    std::pair{1.00, "warn"},
    std::pair{1.00, "error"},
  };

  static constexpr auto overrepresented_cutoff = 0.01;
  // only count first 100k unique sequences
  static constexpr auto max_unique{100'000};
  std::unordered_map<falco_word, std::uint64_t> dups;
  bool add_unique_seqs{true};
  std::uint64_t n_unique{};
  std::uint64_t limit_count{};

  auto
  summarize();

  [[nodiscard]] auto
  format_overrepresented(const std::uint64_t n_reads) const -> std::string;

  [[nodiscard]] auto
  format_duplication_levels() const -> std::string;

  auto
  operator+=(const duplication_results &rhs) -> const duplication_results &;
};

template <>
struct std::formatter<duplication_results> : std::formatter<std::string> {
  auto
  format(const duplication_results &dr, auto &ctx) const {
    return std::formatter<std::string>::format(dr.format_duplication_levels(),
                                               ctx);
  }
};

inline auto
count_seqs_add(const auto seq_itr, const auto sz, duplication_results &dr) {
  const auto fw = falco_word(seq_itr, sz);
  const auto itr = dr.dups.find(fw);
  if (itr != std::cend(dr.dups)) {
    ++itr->second;
    ++dr.limit_count;
    return;
  }
  dr.dups.emplace(fw, 1);
  ++dr.n_unique;
  ++dr.limit_count;
  dr.add_unique_seqs = (dr.n_unique < duplication_results::max_unique);
}

inline auto
count_seqs_dups(const auto seq_itr, const auto sz, duplication_results &dr) {
  const auto fw = falco_word(seq_itr, sz);
  const auto itr = dr.dups.find(fw);
  if (itr != std::cend(dr.dups))
    ++itr->second;
}

inline auto
count_seqs(const auto seq_itr, const auto sz, duplication_results &dr) {
  if (dr.add_unique_seqs)
    count_seqs_add(seq_itr, sz, dr);
  else
    count_seqs_dups(seq_itr, sz, dr);
}

#endif  // SRC_DUPLICATION_RESULTS_HPP_
