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

#ifndef SRC_QUALITY_SCORE_HPP_
#define SRC_QUALITY_SCORE_HPP_

#include <algorithm>
#include <array>
#include <format>
#include <iterator>
#include <ranges>
#include <stdexcept>

namespace falco {
static constexpr auto max_qual_val = 126;
// clang-format off
static constexpr auto min_qual_cutoffs = std::array{
  33,  // Sanger / Illumina 1.9
  64,  // Solexa / Illumina 1.3 / Illumina 1.5
};
static constexpr auto min_qual_offsets = std::array{
  0,
  33,  // Sanger / Illumina 1.9
  64,  // Solexa / Illumina 1.3 / Illumina 1.5
};
static constexpr auto format_labels = std::array{
  "Unknown",
  "Sanger / Illumina 1.9",
  "Solexa / Illumina 1.3 / Illumina 1.5 / Illumina 1.8",
};
// clang-format on
}  // namespace falco

[[nodiscard]] inline auto
identify_quality_score_idx(const auto &qual_counts) {
  // get minimum quality score seen from an array of counts of observed scores
  const auto gt0 = [](const auto c) { return c > 0; };
  const auto min_qual_itr = std::ranges::find_if(qual_counts, gt0);
  auto min_qual = std::distance(std::cbegin(qual_counts), min_qual_itr);
  if (min_qual > falco::max_qual_val)
    throw std::runtime_error(
      std::format("found invalid quality score: {}", min_qual));
  // now find the index among minimum values for formats
  const auto itr = std::ranges::upper_bound(falco::min_qual_cutoffs, min_qual);
  auto idx = std::distance(std::cbegin(falco::min_qual_cutoffs), itr);
  if (idx == std::size(falco::min_qual_cutoffs))
    --idx;
  return idx;
}

[[nodiscard]] inline auto
identify_quality_score_encoding(const auto &qual_counts) {
  return falco::format_labels[identify_quality_score_idx(qual_counts)];
}

[[nodiscard]] inline auto
identify_quality_score_offset(const auto &qual_counts) {
  return falco::min_qual_offsets[identify_quality_score_idx(qual_counts)];
}

#endif  // SRC_QUALITY_SCORE_HPP_
