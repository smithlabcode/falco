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

#include "nlohmann/json.hpp"

#include <algorithm>
#include <array>
#include <format>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

namespace falco {
enum class encoding : std::uint8_t {
  unknown = 0,
  sanger = 1,
  solexa = 2,
};

NLOHMANN_JSON_SERIALIZE_ENUM(  //
  encoding,                    //
  {
    {encoding::unknown, "Unknown"},
    {encoding::sanger, "Sanger / Illumina 1.9"},
    {encoding::solexa, "Solexa / Illumina <= 1.8"},
  })

static constexpr auto bam_qual_offset = 33;
static constexpr auto max_qual_val = 126;
static constexpr auto sanger_min_qual = 33;
static constexpr auto solexa_min_qual = 64;
// clang-format off
static constexpr auto min_qual_offsets = std::array{
  0,
  sanger_min_qual,  // Sanger / Illumina 1.9
  solexa_min_qual,  // Solexa / Illumina 1.3 / Illumina 1.5
};
static constexpr auto format_labels = std::array{
  "Unknown",
  "Sanger / Illumina 1.9",
  "Solexa / Illumina <= 1.8",
};
using qual_array = std::array<std::uint64_t, max_qual_val>;
// clang-format on
}  // namespace falco

[[nodiscard]] static inline auto
set_quality_score_encoding_impl(const auto &qual_counts,
                                const std::int32_t qual_offset) {
  const auto gt0 = [](const auto c) { return c > 0; };
  const auto first_non_zero = [&](const auto &v) {
    return std::distance(std::cbegin(v), std::ranges::find_if(v, gt0));
  };
  const auto min_qual =
    std::ranges::min(qual_counts | std::views::transform(first_non_zero)) +
    qual_offset;
  if (min_qual > falco::max_qual_val)
    throw std::runtime_error("invalid qual score: " + std::to_string(min_qual));
  if (min_qual < falco::sanger_min_qual)
    return falco::encoding::unknown;
  if (min_qual < falco::solexa_min_qual)
    return falco::encoding::sanger;
  return falco::encoding::solexa;
}

[[nodiscard]] inline auto
set_quality_score_encoding(const auto &qual_counts, auto &info) {
  const auto shift = is_mapped_reads(info.format) ? falco::bam_qual_offset : 0;
  info.encoding = set_quality_score_encoding_impl(qual_counts, shift);
}

[[nodiscard]] inline auto
get_quality_score_offset(const auto encoding) {
  // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-constant-array-index)
  return falco::min_qual_offsets[std::to_underlying(encoding)];
}

[[nodiscard]] inline auto
get_quality_score_label(const auto encoding) {
  return falco::format_labels[std::to_underlying(encoding)];
}

inline auto
adjust_fastq_qual_encoding(std::vector<falco::qual_array> &qual_by_pos,
                           falco::qual_array &qual_by_read,
                           const falco::encoding enc) {
  const auto qual_offset = get_quality_score_offset(enc);
  // cppcheck-suppress constParameterReference
  const auto shift_and_fill = [&](auto &x) {
    const auto itr = std::shift_left(std::begin(x), std::end(x), qual_offset);
    std::ranges::fill_n(itr, qual_offset, 0);
  };
  std::ranges::for_each(qual_by_pos, shift_and_fill);
  shift_and_fill(qual_by_read);
}

#endif  // SRC_QUALITY_SCORE_HPP_
