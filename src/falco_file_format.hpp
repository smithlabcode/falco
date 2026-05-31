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

#ifndef SRC_FALCO_FILE_FORMAT_HPP_
#define SRC_FALCO_FILE_FORMAT_HPP_

#include "nlohmann/json.hpp"

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <format>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <memory>  // IWYU pragma: keep
#include <numeric>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>  // IWYU pragma: keep
#include <type_traits>
#include <utility>
#include <vector>

namespace falco {
enum class file_format : std::uint8_t {
  unknown,
  fastq,
  fastq_gz,
  fastq_bgzf,
  sam,
  bam,
};

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

[[nodiscard]] inline auto
get_file_format(const std::string &filename)
  -> std::tuple<falco::file_format, std::string> {
  std::unique_ptr<htsFile, int (*)(htsFile *)> fp(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!fp)
    throw std::runtime_error("failed to open file: " + filename);
  const auto fmt = hts_get_format(fp.get());
  const auto descr = std::string{
    std::unique_ptr<char, void (*)(void *)>(hts_format_description(fmt), &free)
      .get()};
  if (fmt->format == bam)
    return {falco::file_format::bam, descr};
  if (fmt->format == sam)
    return {falco::file_format::sam, descr};
  if (fmt->format == fastq_format) {
    // check for compression
    if (fmt->compression == bgzf)
      return {falco::file_format::fastq_bgzf, descr};
    if (fmt->compression == gzip)
      return {falco::file_format::fastq_gz, descr};
    return {falco::file_format::fastq, descr};
  }
  return {falco::file_format::unknown, descr};
}

#endif  // SRC_FALCO_FILE_FORMAT_HPP_
