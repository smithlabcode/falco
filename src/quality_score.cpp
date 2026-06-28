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

#include "quality_score.hpp"

#include "falco_file_format.hpp"
#include "falco_utils.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <format>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <utility>
#include <vector>

[[nodiscard]] static inline auto
identify_encoding_impl(const auto &qual_counts,
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

[[nodiscard]] auto
identify_encoding(const std::vector<falco::qual_array> &qual_counts,
                  file_info &info) -> falco::encoding {
  const auto shift = is_mapped_reads(info.format) ? falco::bam_qual_offset : 0;
  return identify_encoding_impl(qual_counts, shift);
}

[[nodiscard]] auto
get_quality_score_offset(const falco::encoding e) -> std::int64_t {
  // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-constant-array-index)
  return falco::min_qual_offsets[std::to_underlying(e)];
}

[[nodiscard]] auto
get_quality_score_label(const falco::encoding e) -> std::string {
  // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-constant-array-index)
  return falco::format_labels[std::to_underlying(e)];
}

auto
adjust_fastq_qual_encoding(std::vector<falco::qual_array> &qual_by_pos,
                           falco::qual_array &qual_by_read,
                           const falco::encoding enc) -> void {
  const auto qual_offset = get_quality_score_offset(enc);
  // cppcheck-suppress constParameterReference
  const auto shift_and_fill = [&](auto &x) {
    const auto itr = std::shift_left(std::begin(x), std::end(x), qual_offset);
    std::ranges::fill_n(itr, qual_offset, 0);
  };
  std::ranges::for_each(qual_by_pos, shift_and_fill);
  shift_and_fill(qual_by_read);
}

[[nodiscard]] auto
to_string(const falco::encoding e) -> std::string {
  const auto u = std::to_underlying(e);
  assert(u < std::size(falco::format_labels));
  return falco::format_labels[u];  // NOLINT(*-pro-bounds-constant-array-index)
}
