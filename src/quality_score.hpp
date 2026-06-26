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

#include <array>
#include <cstdint>
#include <format>
#include <string>

struct file_info;

namespace falco {
enum class encoding : std::uint8_t {
  unknown = 0,
  sanger = 1,
  solexa = 2,
};

// NOLINTNEXTLINE
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

[[nodiscard]] auto
to_string(const falco::encoding e) -> std::string;

[[nodiscard]] auto
get_quality_score_offset(const falco::encoding e) -> std::int64_t;

[[nodiscard]] auto
identify_encoding(const std::vector<falco::qual_array> &qual_counts,
                  file_info &info) -> falco::encoding;

[[nodiscard]] auto
get_quality_score_label(const falco::encoding e) -> std::string;

auto
adjust_fastq_qual_encoding(std::vector<falco::qual_array> &qual_by_pos,
                           falco::qual_array &qual_by_read,
                           const falco::encoding enc) -> void;

#endif  // SRC_QUALITY_SCORE_HPP_
