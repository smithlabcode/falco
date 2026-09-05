// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "sam_stdin.hpp"
#include "samrec.hpp"
#include "task_queue.hpp"

#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <iterator>
#include <memory>
#include <ranges>
#include <string>
#include <system_error>
#include <tuple>  // IWYU pragma: keep

sam_stdin::sam_stdin(const std::int64_t buf_size) :
  buffer(buf_size + min_buf_size), cursor{std::begin(buffer)},
  last{std::begin(buffer)} {
  skip_header();
}

[[nodiscard]] auto
estimate_n_reads_sam_stdin(const std::string &)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  static constexpr auto assumed_n_reads = 10'000'000;
  static constexpr auto assumed_read_len = 150;
  static constexpr auto assumed_filesize = 0;
  return {assumed_n_reads, assumed_read_len, assumed_filesize};
}

auto
sam_stdin::skip_header() -> void {
  bool pre_header = true;
  while (cursor == std::cbegin(buffer)) {
    load_next();
    if (pre_header && buffer[0] != '@')
      break;
    pre_header = false;
    const auto lim = std::distance(std::begin(buffer), last);
    for (auto i = 1; i < lim && cursor == std::begin(buffer); ++i)
      if (buffer[i - 1] == '\n' && buffer[i] != '@')
        cursor = std::begin(buffer) + i;
    if (cursor == std::cbegin(buffer)) {
      // corner case of '\n' and non-@ in different buffers
      last = std::end(buffer);
      cursor = std::prev(last);
      shift_output_buffer();  // in prep for subsequent load-next
    }
  }
}

auto
sam_stdin::get_chunks(const std::int64_t n_chunks,  //
                      const std::int32_t file_id,   //
                      task_queue &tq,               //
                      std::atomic_int32_t &n_tasks) -> void {
  assert(n_chunks > 0);
  const auto beg_itr = std::begin(buffer);
  const auto end_itr = last;
  // clang-format off
  const auto fwd_to_read_start = [&](auto p) {
    if (p == beg_itr) return p;
    while (p != end_itr && *p != '\n') ++p;
    if (p != end_itr) ++p;
    return p;
  };
  const auto rev_to_read_start = [&](auto p) {
    while (p != beg_itr && *(p - 1) != '\n') --p;
    return p;
  };
  // clang-format on
  const std::int64_t n_bytes_available = std::distance(beg_itr, end_itr);
  const auto [chunk_size, remainder] = std::div(n_bytes_available, n_chunks);
  assert(n_chunks > 0);
  auto start_itr = beg_itr;
  auto chunk_end = start_itr;
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(start_itr);
    const auto stop_itr = start_itr + chunk_size + (chunk_idx < remainder);
    chunk_end = fwd_to_read_start(stop_itr);
    if (chunk_idx + 1 == n_chunks) {  // final chunk has only full records
      const auto prev_start = rev_to_read_start(chunk_end);
      const auto n_trailing = std::count(prev_start, end_itr, '\n');
      if (n_trailing == 0)
        chunk_end = prev_start;
    }
    ++n_tasks;
    tq.push(file_id, sam_task_t(chunk_beg, chunk_end));
    start_itr = stop_itr;
  }
  cursor = chunk_end;
}

auto
sam_stdin::shift_output_buffer() -> void {
  if (cursor == std::cbegin(buffer))  // shifting here does nothing
    return;
  last = std::copy(cursor, last, std::begin(buffer));
  cursor = std::begin(buffer);
}

auto
sam_stdin::load_next() -> void {
  auto space = std::distance(last, std::end(buffer));
  std::int64_t n{1};
  while (space > 0 && (n = read(0, std::to_address(last), space)) != 0) {
    if (n == -1)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "error reading fastq from stdin");
    space -= n;
    last += n;
  }
  hit_eof = (n == 0);
}

auto
sam_stdin::make_tasks(const std::int64_t n_chunks,  //
                      const std::int32_t file_id,   //
                      task_queue &tq,               //
                      std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes more tasks
  shift_output_buffer();
  load_next();
  get_chunks(n_chunks, file_id, tq, n_tasks);
}
