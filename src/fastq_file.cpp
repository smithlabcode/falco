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

#include "fastq_file.hpp"
#include "fastq_buffer.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

[[nodiscard]] auto
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

static auto
cleanup_mmap_fastq(fastq_buffer &fq) {
  if (fq.data == nullptr)
    throw std::runtime_error("failed to cleanup memory map");
  const int rc = munmap(static_cast<void *>(fq.data), fq.sz);
  if (rc)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "cleanup_mmap_fastq");
  fq = {nullptr, 0};
}

fastq_file::fastq_file(const std::string &fastq_filename,
                       const std::int64_t buf_size) :
  buf_size{buf_size}, filesize{static_cast<std::int64_t>(
                        std::filesystem::file_size(fastq_filename))} {
  // ADS: file descriptor will be reused
  fd = open(std::data(fastq_filename), O_RDONLY, 0);
  if (fd < 0)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            R"(failed to open file: )" + fastq_filename);
  stop_offset = buf_size;
}

fastq_file::~fastq_file() {
  if (fq.sz > 0)
    cleanup_mmap_fastq(fq);
  if (fd)
    close(fd);  // done with file descriptor
}

[[nodiscard]] auto
fastq_file::get_next() -> const fastq_buffer & {
  start_offset += cursor;
  cursor = start_offset & mv.offset_mask;
  start_offset = start_offset & mv.page_mask;
  stop_offset = std::min(filesize, start_offset + buf_size);
  if (fq.sz > 0)
    cleanup_mmap_fastq(fq);
  fq = mmap_fastq(fd, start_offset, stop_offset);
  return fq;
}

[[nodiscard]] auto
get_chunks(const fastq_buffer &fq, const std::uint64_t start_idx,
           const std::uint64_t stop_idx, const std::uint64_t n_chunks)
  -> std::vector<std::pair<std::uint64_t, std::uint64_t>> {
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
  std::vector<std::pair<std::uint64_t, std::uint64_t>> chunks(n_chunks);
  std::uint64_t block_beg{};
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
