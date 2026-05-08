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

#ifndef SRC_FASTQ_FILE_FORMAT_HPP_
#define SRC_FASTQ_FILE_FORMAT_HPP_

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
#include <numeric>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

enum class file_format : std::uint8_t {
  unknown,
  fastq,
  fastq_gz,
  bam,
};

template <> struct std::formatter<file_format> : std::formatter<std::string> {
  auto
  format(const file_format &f, auto &ctx) const {
    return std::formatter<std::string>::format(
      std::to_string(std::to_underlying(f)), ctx);
  }
};

[[nodiscard]] consteval auto
is_sequence_data(const auto hts_fp) -> bool {
  return htsFormatCategory(hts_get_format(hts_fp)->format) == sequence_data;
}

[[nodiscard]] consteval auto
is_sequence_data(const std::string &filename) -> bool {
  std::unique_ptr<htsFile, int (*)(htsFile *)> fp(
    hts_open(std::data(filename), "r"), &hts_close);
  return is_sequence_data(fp.get());
}

[[nodiscard]] inline auto
get_file_format(const std::string &filename)
  -> std::tuple<file_format, std::string> {
  std::unique_ptr<htsFile, int (*)(htsFile *)> fp(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!fp)
    throw std::runtime_error("failed to open file: " + filename);
  const auto fmt = hts_get_format(fp.get());
  const auto descr = std::string{
    std::unique_ptr<char, void (*)(void *)>(hts_format_description(fmt), &free)
      .get()};
  // check for BAM/SAM
  if (fmt->format == bam || fmt->format == sam)
    return {file_format::bam, descr};
  // check for FASTQ
  if (fmt->format == fastq_format) {
    // check for compressed
    if (fmt->compression == gzip || fmt->compression == bgzf)
      return {file_format::fastq_gz, descr};
    return {file_format::fastq, descr};
  }
  return {file_format::unknown, descr};
}

#endif  // SRC_FASTQ_FILE_FORMAT_HPP_
