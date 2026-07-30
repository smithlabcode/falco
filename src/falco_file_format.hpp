// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_FILE_FORMAT_HPP_
#define SRC_FALCO_FILE_FORMAT_HPP_

#include "nlohmann/json.hpp"

#include <htslib/sam.h>

#include <cstdint>
#include <format>
#include <iterator>
#include <memory>
#include <string>
#include <thread>
#include <tuple>
#include <utility>
#include <variant>

namespace falco {
enum class file_format : std::uint8_t {
  unknown,
  fastq,
  fastq_gz,
  fastq_bgzf,
  sam,
  bam,
};

// NOLINTNEXTLINE
NLOHMANN_JSON_SERIALIZE_ENUM(  //
  file_format,                 //
  {
    {file_format::unknown, "unknown"},
    {file_format::fastq, "fastq"},
    {file_format::fastq_gz, "fastq_gz"},
    {file_format::fastq_bgzf, "fastq_bgzf"},
    {file_format::sam, "SAM"},
    {file_format::bam, "BAM"},
  })

[[nodiscard]] constexpr inline auto
is_mapped_reads(const file_format f) {
  return f == file_format::sam || f == file_format::bam;
}

[[nodiscard]] constexpr inline auto
is_bam(const file_format f) {
  return f == file_format::bam;
}

[[nodiscard]] constexpr inline auto
is_bgzf(const file_format f) {
  return f == file_format::bam || f == file_format::fastq_bgzf;
}

[[nodiscard]] constexpr inline auto
is_plain(const file_format f) {
  return f == file_format::sam || f == file_format::fastq;
}

}  // namespace falco

template <>
struct std::formatter<falco::file_format> : std::formatter<std::string> {
  auto
  format(const falco::file_format &f, auto &ctx) const {
    return std::formatter<std::string>::format(
      std::to_string(std::to_underlying(f)), ctx);
  }
};

// ADS: unused?
[[nodiscard]] constexpr auto
is_sequence_data(const auto hts_fp) -> bool {
  return hts_get_format(hts_fp)->category == sequence_data;
}

[[nodiscard]] constexpr auto
is_sequence_data(const std::string &filename) -> bool {
  using hts_file_unique_ptr = std::unique_ptr<htsFile, int (*)(htsFile *)>;
  hts_file_unique_ptr fp(hts_open(std::data(filename), "r"), &hts_close);
  return is_sequence_data(fp.get());
}

[[nodiscard]] auto
get_file_format(const std::string &filename)
  -> std::tuple<falco::file_format, std::string>;

[[nodiscard]] auto
remove_extension(const std::string &filename) -> std::string;

#endif  // SRC_FALCO_FILE_FORMAT_HPP_
