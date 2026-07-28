// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "fastq_gz_file.hpp"
#include "falco_utils.hpp"
#include "task_queue.hpp"

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

[[nodiscard]] auto
estimate_n_reads_fastq_gz(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  static constexpr auto fastq_lines_per_read = 4;
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
  const auto read_len_est =
    estimate_read_length_fastq_chunk(buf, std::size(buf));
  return {static_cast<std::uint64_t>(n_reads_est), read_len_est, filesize};
}

auto
fastq_gz_file::get_chunks(const std::int64_t n_chunks,  //
                          const std::int32_t file_id,   //
                          task_queue &tq,               //
                          std::atomic_int32_t &n_tasks  //
                          ) -> void {
  static constexpr auto rec_lines = 4;  // FASTQ
  assert(n_chunks > 0);
  const auto data = std::data(outbuf);
  const auto not_read_start = [](const auto s, const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return s[p] != '@' || (p > 0 && s[p - 1] != '\n') ||
           (p > 2 && s[p - 2] == '+' && s[p - 3] == '\n');
  };
  const auto fwd_to_read_start = [&](const auto s, auto p) {
    if (p == 0)
      return p;
    while (p < buf_sz && not_read_start(s, p))
      ++p;
    return p;
  };
  const auto rev_to_read_start = [&](const auto s, auto p) {
    while (p > 0 && (p == buf_sz || not_read_start(s, p)))
      --p;
    return p;
  };
  const auto n_bytes_available = buf_sz - cursor;
  const auto [chunk_size, remainder] = std::div(n_bytes_available, n_chunks);
  assert(n_chunks > 0);
  std::int64_t start_pos = cursor;
  std::int64_t chunk_end{};
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(buf_data, start_pos);
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    chunk_end = fwd_to_read_start(buf_data, stop_pos);
    if (chunk_idx == n_chunks - 1) {
      // make sure final chunk includes only full records
      const auto prev_start = rev_to_read_start(buf_data, chunk_end);
      const auto trailing_lines = std::count(
        std::cbegin(outbuf) + prev_start, std::cbegin(outbuf) + buf_sz, '\n');
      if (trailing_lines < rec_lines)
        chunk_end = prev_start;
    }
    ++n_tasks;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    tq.push(file_id, fq_task_t(data + chunk_beg, data + chunk_end));
    start_pos = stop_pos;
  }
  cursor = chunk_end;
}

auto
make_tasks(fastq_gz_file &reads_file,    //
           const std::int64_t n_chunks,  //
           const std::int32_t file_id,   //
           task_queue &tq,               //
           std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes tasks
  reads_file.load_next();
  reads_file.get_chunks(n_chunks, file_id, tq, n_tasks);
}
