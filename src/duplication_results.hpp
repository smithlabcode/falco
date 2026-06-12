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

#include <array>
#include <cstdint>
#include <string>
#include <utility>

struct duplication_results {
  static constexpr auto max_n_reads_total{1'000'000};
  static constexpr auto default_read_skip{10};
  static constexpr auto overrep_cutoff = 0.001;
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.50, "fail"},
    std::pair{0.70, "warn"},
    std::pair{1.00, "pass"},
  };
  static constexpr auto overrep_grade_cutoffs = std::array{
    std::pair{0.10, "pass"},
    std::pair{1.00, "warn"},
    std::pair{1.00, "fail"},
  };
  std::int32_t read_skip{default_read_skip};
  std::int32_t read_idx{};
  boost::unordered_flat_map<falco_word, std::uint64_t> dups;

  auto
  initialize(const std::uint64_t est_n_reads) -> void;

  [[nodiscard]] auto
  format_overrepresented(std::string &grade) const -> std::string;

  [[nodiscard]] auto
  format_duplication_levels(std::string &grade) const -> std::string;

  auto
  operator+=(const duplication_results &rhs) -> const duplication_results &;

  auto
  count_seqs(const auto seq_itr, const auto sz) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      ++dups[falco_word(seq_itr, sz)];
    }
  }
};

#endif  // SRC_DUPLICATION_RESULTS_HPP_
