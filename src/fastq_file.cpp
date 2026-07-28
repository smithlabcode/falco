// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "fastq_file.hpp"
#include "falco_utils.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdint>
#include <string>
#include <system_error>

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
    return {0LU, 0LU, filesize};

  const auto part_size =
    filesize < n_parts * max_part_size ? filesize / n_parts : max_part_size;

  auto total_newlines = 0ul;
  auto read_len_est = 0ul;
  for (auto i = 0; i < n_parts; ++i) {
    const auto offset = (i * part_size) & page_mask;
    char *data = static_cast<char *>(
      mmap(nullptr, part_size, PROT_READ, MAP_PRIVATE, fd, offset));
    if (data == MAP_FAILED)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to mmap file");
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    total_newlines += std::ranges::count(data, data + part_size, '\n');
    read_len_est += estimate_read_length_fastq_chunk(data, part_size);

    if (munmap(static_cast<void *>(data), part_size))
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to unmap memory");
  }
  close(fd);

  read_len_est /= n_parts;
  const auto n_reads_est =
    as_frac(total_newlines, fastq_lines_per_read) *
    as_frac(static_cast<double>(filesize), (part_size * n_parts));

  return {static_cast<std::uint64_t>(n_reads_est), read_len_est, filesize};
}
