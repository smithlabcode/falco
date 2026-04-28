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

#include <unistd.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <cstdint>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <system_error>

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

template <typename T>
[[nodiscard]] static inline auto
get_chunks(const fastq_buffer &fq, const T start_idx, const T stop_idx,
           const T n_chunks) -> std::vector<std::pair<T, T>> {
  const auto to_read_start = [](const fastq_buffer &fq, auto x) {
    const auto not_read_start = [](const auto &s, const auto p) {
      return s[p] != '@' || s[p - 1] != '\n' ||
             (s[p - 2] == '+' && s[p - 3] == '\n');
    };
    while (x < fq.sz && not_read_start(fq.data, x))
      ++x;
    return x;
  };
  const auto n_elements = stop_idx - start_idx;
  const auto q = n_elements / n_chunks;
  const auto r = n_elements - q * n_chunks;
  std::vector<std::pair<T, T>> chunks(n_chunks);
  T block_beg{};
  for (auto i = 0u; i < n_chunks; ++i) {
    auto chunk_beg = start_idx + block_beg;
    chunk_beg = chunk_beg == 0 ? chunk_beg : to_read_start(fq, chunk_beg);
    const auto sz = i < r ? q + 1 : q;
    const auto block_end = block_beg + sz;
    const auto chunk_end =
      to_read_start(fq, std::max(chunk_beg, start_idx + block_end));
    chunks[i] = {chunk_beg, chunk_end};
    block_beg = block_end;
  }
  return chunks;
}

#endif  // SRC_FASTQ_FILE_HPP_
