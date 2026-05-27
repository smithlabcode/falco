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

#ifndef SRC_CONTAMINANTS_HPP_
#define SRC_CONTAMINANTS_HPP_

#include <algorithm>
#include <iterator>
#include <ranges>
#include <string>
#include <utility>  // IWYU pragma: keep
#include <vector>

extern std::vector<std::pair<std::string, std::string>> contaminants;

void
load_contaminants(const std::string &filename);

// get the longest substring of left that is a prefix of right
[[nodiscard]] static inline auto
get_overlap(const auto &left, const auto &right) {
  const auto left_beg = std::cbegin(left);
  const auto left_end = std::cend(left);
  auto best_n_matches = 0l;
  for (auto left_itr = left_beg; left_itr != left_end; ++left_itr) {
    const auto [mm_left, _] =
      std::ranges::mismatch(std::ranges::subrange(left_itr, left_end), right);
    const auto n_matches = std::distance(left_itr, mm_left);
    best_n_matches = std::max(best_n_matches, n_matches);
  }
  return best_n_matches;
}

[[nodiscard]] inline auto
match_contaminant(const auto &query, const auto &contams) -> std::string {
  static constexpr auto no_hit_label = "No Hit";
  std::string best_name;
  auto best_match = 0l;
  auto best_match_len = 0ul;
  for (const auto &[name, seq] : contams) {
    const auto n_match =
      std::max(get_overlap(query, seq), get_overlap(seq, query));
    if (n_match > best_match) {
      best_name = name;
      best_match = n_match;
      best_match_len = std::size(seq);
    }
  }
  const auto match_cutoff = std::min(best_match_len, std::size(query)) / 2.0;
  // If any sequence is a match, return the best one
  return best_match < match_cutoff ? no_hit_label : best_name;
}

#endif  // SRC_CONTAMINANTS_HPP_
