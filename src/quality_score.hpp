// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_QUALITY_SCORE_HPP_
#define SRC_QUALITY_SCORE_HPP_

#include "nlohmann/json.hpp"

#include <array>
#include <cstdint>
#include <format>
#include <string>
#include <vector>  // IWYU pragma: keep

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
using qual_array = std::array<std::uint64_t, max_qual_val + 1>;
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

[[nodiscard]] static inline auto
count_quals(auto qual_itr, const auto qual_end,
            auto &tab) {  // cppcheck-suppress constParameterReference
  auto out_itr = std::begin(tab);
  auto qual_tot = 0;
  while (qual_itr != qual_end) {
    const auto q = *qual_itr++;
    assert(q >= 0 && q <= falco::max_qual_val);
    ++(*out_itr++)[q];
    qual_tot += q;
  }
  return qual_tot;
}

[[nodiscard]] static inline auto
count_quals_rev(auto qual_itr, const auto qual_end,
                auto &tab) {  // cppcheck-suppress constParameterReference
  auto out_itr = std::begin(tab) + std::distance(qual_itr, qual_end);
  auto qual_tot = 0;
  while (qual_itr != qual_end) {
    const auto q = *qual_itr++;
    assert(q >= 0 && q <= falco::max_qual_val);
    ++(*(--out_itr))[q];
    qual_tot += q;
  }
  return qual_tot;
}

static inline auto
count_quals_itr(auto qual_itr, const auto qual_end, auto tab_itr) {
  while (qual_itr != qual_end) {
    assert(*qual_itr >= 0 && *qual_itr <= falco::max_qual_val);
    tab_itr->first += *qual_itr++;
    ++tab_itr->second;
    ++tab_itr;
  }
}

static inline auto
count_quals_itr_rev(auto qual_itr, const auto qual_end, auto tab_itr) {
  tab_itr += std::distance(qual_itr, qual_end);
  while (qual_itr != qual_end) {
    assert(*qual_itr >= 0 && *qual_itr <= falco::max_qual_val);
    --tab_itr;
    tab_itr->first += *qual_itr++;
    ++tab_itr->second;
  }
}

#endif  // SRC_QUALITY_SCORE_HPP_
