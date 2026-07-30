// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "falco_file_format.hpp"

#include <array>
#include <cstdlib>
#include <filesystem>
#include <stdexcept>

[[nodiscard]] auto
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
  if (std::filesystem::file_size(filename) == 0)
    return {falco::file_format::fastq, descr};
  return {falco::file_format::unknown, descr};
}

[[nodiscard]] auto
remove_extension(const std::string &filename) -> std::string {
  using std::literals::string_literals::operator""s;
  // clang-format off
  static const auto extensions = std::array{
    ".fq"s,
    ".fq.gz"s,
    ".fastq"s,
    ".fastq.gz"s,
    ".bam"s,
    ".sam"s,
  };
  // clang-format off
  for (const auto &e : extensions) {
    const std::int64_t non_extn_sz = std::ssize(filename) - std::ssize(e);
    if (filename.ends_with(e) && non_extn_sz > 0)
      return filename.substr(0, non_extn_sz);
  }
  return filename;
}
