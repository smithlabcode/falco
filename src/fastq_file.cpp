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

#include <htslib/hfile.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <memory>
#include <ranges>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename) -> std::uint64_t {
  static constexpr auto n_parts = 10;
  static constexpr auto max_part_size = 1024 * 1024;
  static const auto page_mask = ~(sysconf(_SC_PAGESIZE) - 1);
  int fd{open(std::data(filename), O_RDONLY, 0)};
  if (fd < 0)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
  struct stat buf;
  fstat(fd, &buf);
  const auto filesize = buf.st_size;

  const auto part_size =
    filesize < n_parts * max_part_size ? filesize / n_parts : max_part_size;

  auto total_newlines = 0ul;
  for (auto i = 0; i < n_parts; ++i) {
    const auto offset = (i * part_size) & page_mask;
    char *data = static_cast<char *>(
      mmap(nullptr, part_size, PROT_READ, MAP_PRIVATE, fd, offset));
    if (data == MAP_FAILED)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to mmap file");
    total_newlines += std::ranges::count(data, data + part_size, '\n');
    [[maybe_unused]] const int rc =
      munmap(static_cast<void *>(data), part_size);
  }
  close(fd);

  return (total_newlines / 4.0) *
         (static_cast<double>(filesize) / (part_size * n_parts));
}

[[nodiscard]] auto
estimate_n_reads_fastq_gz(const std::string &filename) -> std::uint64_t {
  static constexpr auto n_bytes = 1024 * 1024;
  std::unique_ptr<BGZF, int (*)(BGZF *)> f(bgzf_open(std::data(filename), "r"),
                                           &bgzf_close);
  std::vector<std::uint8_t> buf(n_bytes);
  const auto r = bgzf_read(f.get(), std::data(buf), n_bytes);
  if (r < 0)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed reading gz input file");
  // ADS: 'htell' function works below because 'f' has no threadpool
  const auto n_compressed_bytes = htell(f.get()->fp);
  const auto total_newlines = std::ranges::count(buf, '\n');
  const double inflation_factor =
    static_cast<double>(n_bytes) / static_cast<double>(n_compressed_bytes);
  const auto filesize = std::filesystem::file_size(filename);
  const auto estimated_uncompressed_file_size =
    std::ceil(inflation_factor * filesize);
  return ((total_newlines / 4.0) / n_bytes) * estimated_uncompressed_file_size;
}
