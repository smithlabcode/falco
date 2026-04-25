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

#ifndef FASTQ_FILE_HPP_
#define FASTQ_FILE_HPP_

#include "fastq_buffer.hpp"

#include <fcntl.h>     // open, O_RDONLY
#include <sys/mman.h>  // mmap, munmap, MAP_FAILED, MAP_PRIVATE
#include <unistd.h>    // close

#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>
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

struct fastq_file {
  std::int64_t buf_size{};
  std::int64_t filesize{};
  fastq_buffer fq{};
  std::int64_t start_offset{};
  std::int64_t stop_offset{};
  std::int64_t cursor{};
  mmap_vals mv{};
  int fd{};

  fastq_file(const std::string &fastq_filename, const std::int64_t buf_size);
  ~fastq_file();

  [[nodiscard]] operator bool() const {
    return stop_offset - start_offset == buf_size;
  }

  [[nodiscard]] auto
  get_next() -> const fastq_buffer &;
};

[[nodiscard]] auto
get_chunks(const fastq_buffer &fq, const std::uint64_t start_idx,
           const std::uint64_t stop_idx, const std::uint64_t n_chunks)
  -> std::vector<std::pair<std::uint64_t, std::uint64_t>>;

#endif  // FASTQ_FILE_HPP_
