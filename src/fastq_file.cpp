// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "fastq_file.hpp"
#include "falco_utils.hpp"
#include "fqrec.hpp"
#include "task_queue.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <compare>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iterator>
#include <memory>
#include <ranges>
#include <span>
#include <string>
#include <system_error>
#include <tuple>
#include <vector>

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  static constexpr auto fastq_lines_per_read = 4;
  static constexpr auto n_parts = 10L;
  static constexpr auto max_part_size = 1024 * 1024;
  std::vector<char> buffer(max_part_size);

  std::unique_ptr<std::FILE, int (*)(std::FILE *)> in(
    std::fopen(std::data(filename), "r"), &std::fclose);
  if (!in)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
  struct stat buf{};
  fstat(fileno(in.get()), &buf);
  const auto filesize = buf.st_size;
  if (filesize < n_parts)
    return {{}, {}, filesize};

  const auto [part_size, remainder] = (filesize < n_parts * max_part_size)
                                        ? std::div(filesize, n_parts)
                                        : std::ldiv_t{max_part_size, 0};
  auto n_lines = 0LU;
  auto readlen_est = 0LU;
  auto offset = 0L;
  for (auto i = 0; i < n_parts; ++i) {
    if (std::fseek(in.get(), offset, SEEK_SET))
      std::system_error(std::make_error_code(std::errc(errno)),
                        "error reading fastq file: " + filename);
    const auto r = std::fread(std::data(buffer), 1, part_size, in.get());
    if (std::ferror(in.get()))
      std::system_error(std::make_error_code(std::errc(errno)),
                        "error reading fastq file: " + filename);
    n_lines += std::ranges::count(std::span{std::cbegin(buffer), r}, '\n');
    readlen_est += estimate_read_length_fastq_chunk(buffer, r);
    offset += part_size + (i < remainder);
  }

  readlen_est /= n_parts;
  const auto n_reads_est =
    as_frac(n_lines, fastq_lines_per_read) *
    as_frac(static_cast<double>(filesize), (part_size * n_parts));

  return {static_cast<std::uint64_t>(n_reads_est), readlen_est, filesize};
}

fastq_file::fastq_file(const std::string &filename,
                       const std::int64_t target_length) :
  target_length{
    std::max(target_length, static_cast<std::int64_t>(min_buf_size))},
  filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
  fd{open(std::data(filename), O_RDONLY, 0)} {
  if (fd < 0)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
}

static inline auto
mmap_fastq(const int fd, const std::int64_t offset, const std::int64_t length,
           auto &data) {
  static constexpr auto prot = PROT_READ;
  static constexpr auto flags = MAP_PRIVATE;
  data = static_cast<char *>(mmap(nullptr, length, prot, flags, fd, offset));
  if (data == MAP_FAILED)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to mmap file");
}

static inline auto
cleanup_mmap_fastq(auto &data, std::int64_t &length) {
  if (data == nullptr)
    return;
  munmap(static_cast<void *>(data), length);
  data = nullptr;
  length = 0;
}

auto
fastq_file::reset() -> void {
  cleanup_mmap_fastq(mmap_data, length);
}

auto
fastq_file::load_next() -> void {
  // memory mapped data is page aligned but the data we need is not
  std::int64_t offset_in_buf = std::distance(mmap_data, std::to_address(last));
  static const auto page_mask = sysconf(_SC_PAGESIZE) - 1;
  std::tie(start_in_file, offset_in_buf) = [&] {
    const auto pos_in_file = start_in_file + offset_in_buf;
    return std::tuple(pos_in_file & (~page_mask), pos_in_file & page_mask);
  }();
  cleanup_mmap_fastq(mmap_data, length);
  stop_in_file = std::min(filesize, start_in_file + target_length);
  if (start_in_file < stop_in_file) {  // this exist for empty files
    length = stop_in_file - start_in_file;
    mmap_fastq(fd, start_in_file, length, mmap_data);
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    buffer = std::span(mmap_data + offset_in_buf, length - offset_in_buf);
  }
}

auto
fastq_file::get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
                       task_queue &tq, std::atomic_int32_t &n_tasks) -> void {
  static constexpr auto rec_lines = 4;  // FASTQ
  assert(n_chunks > 0);
  const auto beg_itr = std::begin(buffer);
  const auto end_itr = std::cend(buffer);
  // clang-format off
  const auto not_read_start = [&](const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return *p != '@' || (p > beg_itr && *(p - 1) != '\n') ||
      (p > beg_itr + 2 && *(p - 2) == '+' && *(p - 3) == '\n');
  };
  const auto forward_to_start = [&](auto pos) {
    if (pos == beg_itr) return pos;
    while (pos < end_itr && not_read_start(pos)) ++pos;
    return pos;
  };
  const auto reverse_to_start = [&](auto pos) {
    while (pos > beg_itr && (pos == end_itr || not_read_start(pos))) --pos;
    return pos;
  };
  // clang-format on
  const auto [chunk_size, remainder] = std::div(std::ssize(buffer), n_chunks);
  auto start_pos = beg_itr;
  auto chunk_end = start_pos;
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = forward_to_start(start_pos);
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    chunk_end = forward_to_start(stop_pos);
    if (chunk_idx + 1 == n_chunks)
      if (const auto prev = reverse_to_start(chunk_end);
          std::count(prev, end_itr, '\n') < rec_lines)
        chunk_end = prev;
    ++n_tasks;
    tq.push(file_id,
            fq_task_t(std::to_address(chunk_beg), std::to_address(chunk_end)));
    start_pos = stop_pos;
  }
  last = chunk_end;
}
