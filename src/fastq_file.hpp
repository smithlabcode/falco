/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SRC_FASTQ_FILE_HPP_
#define SRC_FASTQ_FILE_HPP_

#include "fastq_record.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <cstdint>
#include <filesystem>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

struct mmap_vals {
  // ADS: need to check that buffer size is larger than page size
  const std::int64_t page_size{};
  // ADS: for below, I vaguely recall these values already exist as
  // constants
  const std::int64_t offset_mask{};
  const std::int64_t page_mask{};
  mmap_vals() :
    page_size{sysconf(_SC_PAGESIZE)}, offset_mask{page_size - 1},
    page_mask{~offset_mask} {}
};

[[nodiscard]] inline auto
mmap_fastq(const int fd, const std::int64_t start_offset,
           const std::int64_t stop_offset) -> fastq_buffer {
  const auto n_bytes = stop_offset - start_offset;
  char *data = static_cast<char *>(
    mmap(nullptr, n_bytes, PROT_READ, MAP_PRIVATE, fd, start_offset));
  if (data == MAP_FAILED)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to mmap file");
  return {data, n_bytes};
}

static inline auto
cleanup_mmap_fastq(fastq_buffer &fq) {
  if (fq.data == nullptr)
    return;
  [[maybe_unused]] const int rc = munmap(static_cast<void *>(fq.data), fq.sz);
  fq = {nullptr, 0};
}

struct fastq_file {
  std::int64_t buf_size{};
  std::int64_t filesize{};
  fastq_buffer fq{};
  std::int64_t start_offset{};
  std::int64_t stop_offset{};
  std::int64_t cursor{};
  mmap_vals mv{};
  int fd{};

  fastq_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    fd{open(std::data(filename), O_RDONLY, 0)}, stop_offset{buf_size} {
    if (fd < 0)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              R"(failed to open file: )" + filename);
  }

  ~fastq_file() {
    if (fq.sz > 0)
      cleanup_mmap_fastq(fq);
    if (fd)
      close(fd);  // done with file descriptor
  }

  [[nodiscard]] operator bool() const {
    return stop_offset - start_offset == buf_size;
  }

  [[nodiscard]] auto
  get_next() -> const fastq_buffer & {
    start_offset += cursor;
    cursor = start_offset & mv.offset_mask;
    start_offset = start_offset & mv.page_mask;
    stop_offset = std::min(filesize, start_offset + buf_size);
    if (fq.sz > 0)
      cleanup_mmap_fastq(fq);
    fq = mmap_fastq(fd, start_offset, stop_offset);
    return fq;
  }
};

[[nodiscard]] static inline auto
get_chunks(fastq_file &fq, const fastq_buffer &buf, const std::int64_t n_chunks)
  -> std::vector<std::pair<std::int64_t, std::int64_t>> {
  static constexpr auto lines_per_record = 4;
  const auto not_read_start = [](const auto &s, const auto p) {
    return s[p] != '@' || s[p - 1] != '\n' ||
           (s[p - 2] == '+' && s[p - 3] == '\n');
  };
  const auto forward_to_read_start = [&](const fastq_buffer &buf, auto x) {
    while (x < buf.sz && not_read_start(buf.data, x))
      ++x;
    return x;
  };
  const auto backward_to_read_start = [&](const fastq_buffer &buf, auto x) {
    while (x == buf.sz || not_read_start(buf.data, x))
      --x;
    return x;
  };
  const auto start_idx = fq.cursor;
  const auto stop_idx = buf.sz;
  const auto n_elements = stop_idx - start_idx;
  const auto q = n_elements / n_chunks;
  const auto r = n_elements - q * n_chunks;
  std::vector<std::pair<std::int64_t, std::int64_t>> chunks(n_chunks);
  std::int64_t block_beg{};
  for (auto i = 0u; i < n_chunks; ++i) {
    auto chunk_beg = start_idx + block_beg;
    chunk_beg =
      chunk_beg == 0 ? chunk_beg : forward_to_read_start(buf, chunk_beg);
    const auto sz = i < r ? q + 1 : q;
    const auto block_end = block_beg + sz;
    const auto chunk_end =
      forward_to_read_start(buf, std::max(chunk_beg, start_idx + block_end));
    chunks[i] = {chunk_beg, chunk_end};
    block_beg = block_end;
  }
  const auto prev_read_start =
    backward_to_read_start(buf, chunks.back().second);
  const auto lines_remaining =
    std::ranges::count(buf.data + prev_read_start, buf.data + buf.sz, '\n');
  if (lines_remaining < lines_per_record)
    chunks.back().second = prev_read_start;
  fq.cursor = chunks.back().second;
  return chunks;
}

#endif  // SRC_FASTQ_FILE_HPP_
