// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "fastq_stdin.hpp"
#include "falco_utils.hpp"
#include "fqrec.hpp"
#include "task_queue.hpp"

#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <ranges>
#include <span>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

fastq_stdin::fastq_stdin(const std::int64_t buf_size) :
  buffer(buf_size + min_buf_size), cursor{std::begin(buffer)},
  last{std::begin(buffer)} {}

[[nodiscard]] auto
estimate_n_reads_fastq_stdin(const std::string &)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  static constexpr auto assumed_n_reads = 10'000'000;
  static constexpr auto assumed_read_len = 150;
  static constexpr auto assumed_filesize = 0;
  return {assumed_n_reads, assumed_read_len, assumed_filesize};
}

auto
fastq_stdin::get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
                        task_queue &tq, std::atomic_int32_t &n_tasks) -> void {
  static constexpr auto rec_lines = 4;  // FASTQ
  const auto beg_itr = std::begin(buffer);
  const auto end_itr = last;
  assert(n_chunks > 0);
  // clang-format off
  const auto not_read_start = [&](const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return *p != '@' || (p > beg_itr && *(p - 1) != '\n') ||
      (p > beg_itr + 2 && *(p - 2) == '+' && *(p - 3) == '\n');
  };
  const auto fwd_to_read_start = [&](auto pos) {
    if (pos == beg_itr) return pos;
    while (pos < end_itr && not_read_start(pos)) ++pos;
    return pos;
  };
  const auto rev_to_read_start = [&](auto pos) {
    while (pos > beg_itr && (pos == end_itr || not_read_start(pos))) --pos;
    return pos;
  };
  // clang-format on
  const std::int64_t n_bytes_available = std::distance(beg_itr, last);
  const auto [chunk_size, remainder] = std::div(n_bytes_available, n_chunks);
  std::vector<std::pair<std::int64_t, std::int64_t>> chunks(n_chunks);
  auto start_itr = beg_itr;
  auto chunk_end = start_itr;
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(start_itr);
    const auto stop_itr = start_itr + chunk_size + (chunk_idx < remainder);
    chunk_end = fwd_to_read_start(stop_itr);
    if (chunk_idx + 1 == n_chunks)
      if (const auto prev = rev_to_read_start(chunk_end);
          std::count(prev, std::end(buffer), '\n') < rec_lines)
        chunk_end = prev;
    ++n_tasks;
    tq.push(file_id,
            fq_task_t(std::to_address(chunk_beg), std::to_address(chunk_end)));
    start_itr = stop_itr;
  }
  cursor = chunk_end;
}

auto
fastq_stdin::shift_output_buffer() -> void {
  if (cursor == std::cbegin(buffer))  // shifting here does nothing
    return;
  last = std::copy(cursor, last, std::begin(buffer));
  cursor = std::begin(buffer);
}

auto
fastq_stdin::load_next() -> void {
  const auto n_bytes = std::distance(last, std::end(buffer));
  const auto r = read(0, std::to_address(last), n_bytes);
  if (r == -1)
    std::system_error(std::make_error_code(std::errc(errno)),
                      "error reading fastq from stdin");
  if (r == 0)
    hit_eof = true;
  last += static_cast<std::int64_t>(r);
}

auto
fastq_stdin::make_tasks(const std::int64_t n_chunks,  //
                        const std::int32_t file_id,   //
                        task_queue &tq,               //
                        std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes more tasks
  shift_output_buffer();
  load_next();
  get_chunks(n_chunks, file_id, tq, n_tasks);
}
