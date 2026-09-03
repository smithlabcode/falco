// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "sam_file.hpp"
#include "samrec.hpp"
#include "task_queue.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <string>
#include <system_error>

sam_file::sam_file(const std::string &filename, const std::int64_t buf_size) :
  buffer(buf_size + min_buf_size),
  in(std::fopen(std::data(filename), "r"), &std::fclose) {
  if (!in)
    std::system_error(std::make_error_code(std::errc(errno)),
                      "failed to read file");
  if (!skip_header())
    std::runtime_error("failed to validate SAM header: " + filename);
  cursor = std::begin(buffer);
  last = std::begin(buffer);
}

[[nodiscard]] auto
sam_file::skip_header() -> bool {
  static constexpr auto max_header_size = 256 * 1024 * 1024;
  std::int32_t n_bytes{};
  while (n_bytes < max_header_size) {
    const auto c = std::fgetc(in.get());
    if (c != '@') {  // confirm header line
      const auto r = std::ungetc(c, in.get());
      if (r != c)
        throw std::runtime_error("failed to parse SAM header");
      break;
    }
    ++n_bytes;
    // skip header line
    while (n_bytes++ < max_header_size && std::fgetc(in.get()) != '\n')
      ;
  }
  return n_bytes < max_header_size;
}

auto
sam_file::shift_output_buffer() -> void {
  if (cursor == std::cbegin(buffer))  // shifting here does nothing
    return;
  last = std::copy(cursor, last, std::begin(buffer));
  cursor = std::begin(buffer);
}

auto
sam_file::get_chunks(const std::int64_t n_chunks,  //
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
sam_file::load_next() -> void {
  const auto n_bytes = std::distance(last, std::end(buffer));
  last += static_cast<std::int64_t>(
    std::fread(std::to_address(last), 1, n_bytes, in.get()));
  if (std::ferror(in.get()))
    std::system_error(std::make_error_code(std::errc(errno)),
                      "error reading SAM file");
}

auto
sam_file::make_tasks(const std::int64_t n_chunks,  //
                     const std::int32_t file_id,   //
                     task_queue &tq,               //
                     std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes more tasks
  shift_output_buffer();
  load_next();
  get_chunks(n_chunks, file_id, tq, n_tasks);
}
