// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "fastq_file.hpp"
#include "falco_utils.hpp"
#include "fqrec.hpp"
#include "task_queue.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>  // IWYU pragma: keep
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <ranges>
#include <span>
#include <string>
#include <system_error>
#include <tuple>  // IWYU pragma: keep
#include <type_traits>
#include <utility>
#include <vector>

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  static constexpr auto fastq_lines_per_read = 4;
  static constexpr auto n_parts = 10;
  static constexpr auto max_part_size = 1024 * 1024;
  static const auto page_mask = ~(sysconf(_SC_PAGESIZE) - 1);
  const int fd = open(std::data(filename), O_RDONLY, 0);
  if (fd < 0)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
  struct stat buf{};
  fstat(fd, &buf);
  const auto filesize = buf.st_size;
  if (filesize < n_parts)
    return {{}, {}, filesize};

  const auto part_size =
    filesize < n_parts * max_part_size ? filesize / n_parts : max_part_size;

  auto total_newlines = 0ul;
  auto readlen_est = 0ul;
  for (auto i = 0; i < n_parts; ++i) {
    const auto offset = (i * part_size) & page_mask;
    auto raw = mmap(nullptr, part_size, PROT_READ, MAP_PRIVATE, fd, offset);
    if (raw == MAP_FAILED)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to mmap file");
    const auto data = std::span(static_cast<char *>(raw), part_size);
    total_newlines += std::ranges::count(data, '\n');
    readlen_est += estimate_read_length_fastq_chunk(std::data(data), part_size);
    if (munmap(raw, part_size))
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to unmap memory");
  }
  close(fd);

  readlen_est /= n_parts;
  const auto n_reads_est =
    as_frac(total_newlines, fastq_lines_per_read) *
    as_frac(static_cast<double>(filesize), (part_size * n_parts));

  return {static_cast<std::uint64_t>(n_reads_est), readlen_est, filesize};
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
cleanup_mmap_fastq(auto &buffer, std::int64_t &buffer_size) {
  if (buffer == nullptr)
    return;
  munmap(static_cast<void *>(buffer), buffer_size);
  buffer = nullptr;
  buffer_size = 0;
}

auto
fastq_file::reset() -> void {
  cleanup_mmap_fastq(buffer, buffer_size);
}

auto
fastq_file::load_next() -> void {
  // memory mapped data is page aligned but the data we need is not
  static const auto page_mask = sysconf(_SC_PAGESIZE) - 1;
  std::tie(start_in_file, cursor) = [&] {
    const auto pos_in_file = start_in_file + cursor;
    return std::tuple(pos_in_file & (~page_mask), pos_in_file & page_mask);
  }();
  stop_in_file = std::min(filesize, start_in_file + target_buffer_size);
  if (buffer_size > 0)
    cleanup_mmap_fastq(buffer, buffer_size);
  if (start_in_file < stop_in_file) {  // this exist for empty files
    buffer_size = stop_in_file - start_in_file;
    mmap_fastq(fd, start_in_file, buffer_size, buffer);
  }
}

auto
fastq_file::get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
                       task_queue &tq, std::atomic_int32_t &n_tasks) -> void {
  static constexpr auto rec_lines = 4;  // FASTQ
  assert(n_chunks > 0);
  const auto buf = std::span(buffer, buffer_size);
  // clang-format off
  const auto not_read_start = [&](const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return buf[p] != '@' || (p > 0 && buf[p - 1] != '\n') ||
           (p > 2 && buf[p - 2] == '+' && buf[p - 3] == '\n');
  };
  const auto fwd_to_read_start = [&](auto pos) {
    if (pos == 0) return pos;
    while (pos < buffer_size && not_read_start(pos)) ++pos;
    return pos;
  };
  const auto rev_to_read_start = [&](auto pos) {
    while (pos > 0 && (pos == buffer_size || not_read_start(pos))) --pos;
    return pos;
  };
  // clang-format on
  const auto n_bytes_available = buffer_size - cursor;
  std::int64_t start_pos = cursor;
  const auto [chunk_size, remainder] = std::div(n_bytes_available, n_chunks);
  std::vector<std::pair<std::int64_t, std::int64_t>> chunks(n_chunks);
  std::int64_t chunk_end{start_pos};
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(start_pos);
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    chunk_end = fwd_to_read_start(stop_pos);
    if (chunk_idx + 1 == n_chunks)
      if (const auto prev = rev_to_read_start(chunk_end);
          std::count(std::cbegin(buf) + prev, std::cend(buf), '\n') < rec_lines)
        chunk_end = prev;
    ++n_tasks;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    tq.push(file_id, fq_task_t(buffer + chunk_beg, buffer + chunk_end));
    start_pos = stop_pos;
  }
  cursor = chunk_end;
}
