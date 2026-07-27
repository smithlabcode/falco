// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "fastq_gz_file.hpp"
#include "falco_utils.hpp"

#include <htslib/bgzf.h>
#include <htslib/hfile.h>

#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <ranges>
#include <string>
#include <system_error>
#include <vector>

static constexpr auto fastq_lines_per_read = 4;

[[nodiscard]] static auto
estimate_read_length(const auto &data, const auto n) {
  assert(n >= 1);
  const auto valid = [](const auto c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
  };
  std::vector<std::int64_t> lines;
  for (auto i = 0u; i + 1 < n; ++i)
    if (data[i] == '\n')
      lines.push_back(i + 1);
  if (std::size(lines) < fastq_lines_per_read)
    return 1ul;
  auto total = 0ul;
  for (const auto l : lines | std::views::adjacent<fastq_lines_per_read - 1>)
    if (data[std::get<0>(l)] == '@' && data[std::get<2>(l)] == '+' &&
        valid(data[std::get<1>(l)]))
      // cppcheck-suppress useStlAlgorithm
      total += (std::get<2>(l) - std::get<1>(l)) - 1;
  return total / (std::size(lines) / fastq_lines_per_read);
}

// ADS: currently this is identical to the function for BGZF
[[nodiscard]] auto
estimate_n_reads_fastq_gz(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  static constexpr auto n_bytes = 1024 * 1024;
  std::unique_ptr<BGZF, int (*)(BGZF *)> f(bgzf_open(std::data(filename), "r"),
                                           &bgzf_close);
  std::vector<std::uint8_t> buf(n_bytes);
  const auto r = bgzf_read(f.get(), std::data(buf), n_bytes);
  if (r < 0)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed reading gz input file");
  // ADS: 'htell' function works below because 'f' has no threadpool
  const auto n_compressed_bytes = htell(f.get()->fp);
  const auto total_newlines = std::ranges::count(buf, '\n');

  const auto inflation_factor = as_frac(n_bytes, n_compressed_bytes);
  const auto filesize = std::filesystem::file_size(filename);
  const auto estimated_uncompressed_file_size =
    inflation_factor * static_cast<double>(filesize);
  const auto n_reads_est = as_frac(total_newlines, fastq_lines_per_read) *
                           as_frac(estimated_uncompressed_file_size, n_bytes);
  const auto read_len_est = estimate_read_length(buf, std::size(buf));
  return {static_cast<std::uint64_t>(n_reads_est), read_len_est, filesize};
}

[[nodiscard]] auto
fastq_gz_file::get_chunks(const std::int64_t n_chunks) -> fq_chunks_t {
  assert(n_chunks > 0);
  static constexpr auto rec_lines = 4;  // FASTQ
  const auto not_read_start = [](const auto s, const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return s[p] != '@' || (p > 0 && s[p - 1] != '\n') ||
           (p > 2 && s[p - 2] == '+' && s[p - 3] == '\n');
  };
  const auto fwd_to_read_start = [&](auto pos) {
    if (pos == 0)
      return pos;
    while (pos < buf_sz && not_read_start(buf_data, pos))
      ++pos;
    return pos;
  };
  const auto rev_to_read_start = [&](auto pos) {
    while (pos > 0 && (pos == buf_sz || not_read_start(buf_data, pos)))
      --pos;
    return pos;
  };
  const auto n_bytes_available = buf_sz - cursor;
  const auto [chunk_size, remainder] = std::div(n_bytes_available, n_chunks);
  std::vector<std::pair<std::int64_t, std::int64_t>> chunks(n_chunks);
  std::int64_t start_pos = cursor;
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(start_pos);
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    const auto chunk_end = fwd_to_read_start(stop_pos);
    chunks[chunk_idx] = {chunk_beg, chunk_end};
    start_pos = stop_pos;
  }
  // make sure final chunk includes only full records
  const auto prev_start = rev_to_read_start(chunks.back().second);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  if (std::count(buf_data + prev_start, buf_data + buf_sz, '\n') < rec_lines)
    chunks.back().second = prev_start;
  cursor = chunks.back().second;

  fq_chunks_t tasks;
  tasks.reserve(std::size(chunks));
  std::ranges::for_each(chunks, [&](const auto &chunk) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    tasks.emplace_back(buf_data + chunk.first, buf_data + chunk.second);
  });
  return tasks;
}
