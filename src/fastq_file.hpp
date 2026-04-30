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

#include <htslib/bgzf.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <memory>
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
cleanup_mmap_fastq(fastq_buffer &buf) {
  if (buf.data == nullptr)
    return;
  [[maybe_unused]] const int rc = munmap(static_cast<void *>(buf.data), buf.sz);
  buf = {nullptr, 0};
}

struct fastq_file {
  std::int64_t buf_size{};
  std::int64_t filesize{};
  fastq_buffer buf{};
  std::int64_t start_offset{};
  std::int64_t stop_offset{};
  std::int64_t cursor{};
  mmap_vals mv{};
  int fd{};

  fastq_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    stop_offset{buf_size}, fd{open(std::data(filename), O_RDONLY, 0)} {
    if (fd < 0)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              R"(failed to open file: )" + filename);
  }

  ~fastq_file() {
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
    if (buf.sz > 0)
      cleanup_mmap_fastq(buf);
    buf = mmap_fastq(fd, start_offset, stop_offset);
    return buf;
  }
};

struct fastq_gz_file {
  std::int64_t buf_size{};
  std::int64_t filesize{};
  fastq_buffer buf{};
  std::int64_t start_offset{};
  std::int64_t stop_offset{};
  std::int64_t cursor{};
  std::unique_ptr<BGZF, int (*)(BGZF *)> f;

  // clang-format off
  fastq_gz_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    stop_offset{buf_size},
    f(bgzf_open(std::data(filename), "r"), &bgzf_close)
  {
    if (!f)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              R"(failed to open file: )" + filename);
    buf.data = new char[buf_size];
  }
  // clang-format on

  ~fastq_gz_file() { delete[] buf.data; }

  [[nodiscard]] operator bool() const { return stop_offset == buf_size; }

  [[nodiscard]] auto
  get_next() -> const fastq_buffer & {
    if (cursor > 0) {
      std::memmove(buf.data, buf.data + cursor, buf.sz - cursor);
      cursor = buf.sz - cursor;
    }
    start_offset = cursor;
    const auto n_bytes = buf_size - start_offset;
    const auto r = bgzf_read(f.get(), buf.data + start_offset, n_bytes);
    if (stop_offset < 0) {
      // ADS: cleanup
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed reading gz input file");
    }
    stop_offset = r;
    stop_offset += start_offset;
    buf.sz = stop_offset;
    cursor = 0;  // cursor always moves back to zero if buffer is not mmapped
    return buf;
  }
};

[[nodiscard]] static inline auto
get_chunks(auto &fq, const fastq_buffer &buf, const std::int64_t n_chunks)
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
