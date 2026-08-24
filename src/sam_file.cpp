// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "sam_file.hpp"
#include "samrec.hpp"
#include "task_queue.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <iterator>
#include <ranges>
#include <span>
#include <stdexcept>
#include <system_error>

sam_file::sam_file(const std::string &filename, const std::int64_t buf_size) :
  buffer(buf_size + min_buf_size),
  in(std::fopen(std::data(filename), "r"), &std::fclose) {
  if (!in)
    std::system_error(std::make_error_code(std::errc(errno)),
                      "failed to read file");
  if (!skip_header())
    std::runtime_error("failed to validated SAM file header: " + filename);
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
  if (cursor == std::cbegin(buffer))
    // shifting at cursor == 0 would do nothing
    return;
  const auto n_bytes = std::ranges::distance(cursor, last);
  std::copy(cursor, last, std::begin(buffer));
  cursor = std::begin(buffer);
  last = cursor + n_bytes;
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
  const auto n_bytes_available = std::distance(beg_itr, end_itr);
  const auto chunk_size = (n_bytes_available + n_chunks - 1) / n_chunks;
  assert(n_chunks > 0);
  auto start_pos = beg_itr;
  auto chunk_end = beg_itr;
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(start_pos);
    auto stop_pos = start_pos + chunk_size;
    chunk_end = fwd_to_read_start(stop_pos);
    if (chunk_idx + 1 == n_chunks) {  // final chunk has only full records
      const auto prev_start = rev_to_read_start(chunk_end);
      const auto n_trailing = std::count(prev_start, end_itr, '\n');
      if (n_trailing == 0)
        chunk_end = prev_start;
    }
    ++n_tasks;
    tq.push(file_id, sam_task_t(chunk_beg, chunk_end));
    start_pos = stop_pos;
  }
  cursor = chunk_end;
}

auto
sam_file::make_tasks(const std::int64_t n_chunks,  //
                     const std::int32_t file_id,   //
                     task_queue &tq,               //
                     std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes more tasks
  shift_output_buffer();
  auto inbuf = std::span(last, std::end(buffer));
  const auto r = std::fread(std::data(inbuf), 1, std::size(inbuf), in.get());
  if (std::ferror(in.get()))
    std::system_error(std::make_error_code(std::errc(errno)),
                      "error reading SAM file");
  last += static_cast<std::int64_t>(r);
  get_chunks(n_chunks, file_id, tq, n_tasks);
}
