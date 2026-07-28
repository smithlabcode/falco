// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "sam_file.hpp"
#include "samrec.hpp"      // IWYU pragma: keep
#include "task_queue.hpp"  // IWYU pragma: keep

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <iterator>

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
  char c{};
  while (n_bytes < max_header_size) {
    c = std::fgetc(in.get());
    if (c != '@') {  // confirm header line
      std::ungetc(c, in.get());
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
  const auto n_bytes = last - cursor;
  std::copy_n(std::cbegin(buffer) + cursor, n_bytes, std::begin(buffer));
  last = n_bytes;
  cursor = 0;
}

auto
sam_file::get_chunks(const std::int64_t n_chunks,  //
                     const std::int32_t file_id,   //
                     task_queue &tq,               //
                     std::atomic_int32_t &n_tasks  //
                     ) -> void {
  assert(n_chunks > 0);
  const auto data = std::data(buffer);
  // clang-format off
  const auto fwd_to_read_start = [&](auto p) {
    if (p == 0) return p;
    while (p < last && data[p] != '\n') ++p;
    if (p < last) ++p;
    return p;
  };
  const auto rev_to_read_start = [&](auto p) {
    while (p > 0 && data[p - 1] != '\n') --p;
    return p;
  };
  // clang-format on
  const auto n_bytes_available = last;
  const auto chunk_size = (n_bytes_available + n_chunks - 1) / n_chunks;
  assert(n_chunks > 0);
  std::int64_t start_pos{};
  std::int64_t chunk_end{};
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(start_pos);
    auto stop_pos = start_pos + chunk_size;
    chunk_end = fwd_to_read_start(stop_pos);
    if (chunk_idx == n_chunks - 1) {  // final chunk has only full records
      const auto prev_start = rev_to_read_start(chunk_end);
      const auto n_trailing = std::count(data + prev_start, data + last, '\n');
      if (n_trailing == 0)
        chunk_end = prev_start;
    }
    ++n_tasks;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    tq.push(file_id, sam_task_t(data + chunk_beg, data + chunk_end));
    start_pos = stop_pos;
  }
  cursor = chunk_end;
}

auto
sam_file::make_tasks(const std::int64_t n_chunks,  //
                     const std::int32_t file_id,   //
                     task_queue &tq,               //
                     std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes tasks

  if (cursor > 0)
    shift_output_buffer();
  last +=
    std::fread(std::data(buffer) + last, 1, std::size(buffer) - last, in.get());
  if (std::ferror(in.get()))
    std::system_error(std::make_error_code(std::errc(errno)),
                      "failed to read file");

  get_chunks(n_chunks, file_id, tq, n_tasks);
}
