// SPDX-License-Identifier: MIT; (c) 2026 Andrew D Smith

#include "fastq_file.hpp"
#include "falco_utils.hpp"

#include <htslib/bgzf.h>  // for BGZF
#include <htslib/hfile.h>

#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#include <cstdint>
#include <filesystem>
#include <memory>
#include <ranges>  // IWYU pragma: keep
#include <string>
#include <system_error>
#include <vector>

static constexpr auto fastq_lines_per_read = 4;

[[nodiscard]] static auto
estimate_read_length(const auto &data, const auto n) {
  assert(n >= 1);
  const auto valid = [](const auto c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N';
  };
  std::vector<std::int64_t> lines;
  for (auto i = 0u; i + 1 < n; ++i)
    if (data[i] == '\n')
      lines.push_back(i + 1);
  if (std::size(lines) < fastq_lines_per_read)
    return 1ul;
  auto total = 0ul;
  for (const auto l : lines | std::views::adjacent<fastq_lines_per_read - 1>)
    if (data[std::get<0>(l)] == '@' && data[std::get<2>(l)] == '+' &&
        valid(data[std::get<1>(l)]))
      // cppcheck-suppress useStlAlgorithm
      total += (std::get<2>(l) - std::get<1>(l)) - 1;
  return total / (std::size(lines) / fastq_lines_per_read);
}

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
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
    // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-pointer-arithmetic)
    total_newlines += std::ranges::count(data, data + part_size, '\n');
    read_len_est += estimate_read_length(data, part_size);

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

[[nodiscard]] auto
estimate_n_reads_fastq_bgzf(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
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

  const auto inflation_factor = as_frac(n_bytes, n_compressed_bytes);
  const auto filesize = std::filesystem::file_size(filename);
  const auto estimated_uncompressed_file_size =
    inflation_factor * static_cast<double>(filesize);
  const auto n_reads_est = as_frac(total_newlines, fastq_lines_per_read) *
                           as_frac(estimated_uncompressed_file_size, n_bytes);
  const auto read_len_est = estimate_read_length(buf, std::size(buf));
  return {static_cast<std::uint64_t>(n_reads_est), read_len_est, filesize};
}

[[nodiscard]] auto
estimate_n_reads_fastq_gz(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  return estimate_n_reads_fastq_bgzf(filename);
}

auto
fastq_bgzf_file::load_next() -> std::vector<bgzf_block_t> {
  if (output_cursor > 0) {
    // make room in the output buffer
    std::copy_n(std::cbegin(output_buffer) + output_cursor,
                output_last - output_cursor, std::begin(input_buffer));
    input_last = (output_last - output_cursor);
    input_cursor = 0;  // because it's output
  }
  br.reset();  // to reuse the preallocated buffer
  std::vector<bgzf_block_t> tasks;
  const std::int64_t outbuf_limit =
    std::size(input_buffer) - max_bgzf_block_size;
  // ADS: below, seems a bad pattern to need to keep a task as a class member to
  // hold result of previous calls
  if (!task && !hit_eof)
    task = br.get_decomp_task();
  std::int64_t n_bytes = input_last;
  while (task && n_bytes < outbuf_limit) {
    const auto task_size = task.size;
    // here is where it is confusing: get the next, use the current, swap
    char *out_itr = std::data(input_buffer) + n_bytes;
    n_bytes += task_size;
    bgzf_block_t next_task = br.get_decomp_task();
    task.out_itr = out_itr;
    tasks.emplace_back(std::move(task));
    task = std::move(next_task);
  }
  input_last = n_bytes;
  if (!task && br.has_out())  // seems like it could fail if the output buffer
                              // is exactly full
    hit_eof = true;
  return tasks;
}

[[nodiscard]] auto
fastq_bgzf_file::get_chunks(const std::int64_t n_chunks) -> fq_chunks_t {
  static constexpr auto rec_lines = 4;  // FASTQ
  assert(n_chunks > 0);

  std::swap(input_buffer, output_buffer);
  output_last = input_last;
  output_cursor = input_cursor;

  const auto data = std::data(output_buffer);
  const auto not_read_start = [](const auto s, const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return s[p] != '@' || (p > 0 && s[p - 1] != '\n') ||
           (p > 2 && s[p - 2] == '+' && s[p - 3] == '\n');
  };
  const auto fwd_to_read_start = [&](auto pos) {
    if (pos == 0)
      return pos;
    while (pos < output_last && not_read_start(data, pos))
      ++pos;
    return pos;
  };
  const auto rev_to_read_start = [&](auto pos) {
    while (pos > 0 && (pos == output_last || not_read_start(data, pos)))
      --pos;
    return pos;
  };
  const auto n_bytes_available = output_last - output_cursor;
  const auto [chunk_size, remainder] = std::div(n_bytes_available, n_chunks);
  fq_chunks_t chunks(n_chunks);
  std::int64_t start_pos = output_cursor;
  std::int64_t chunk_end{};
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(start_pos);
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    chunk_end = fwd_to_read_start(stop_pos);
    chunks[chunk_idx] = {data + chunk_beg, data + chunk_end};
    start_pos = stop_pos;
  }
  // make sure final chunk includes only full records
  const auto prev_start = rev_to_read_start(chunk_end);
  const auto trailing_lines =
    std::count(data + prev_start, data + output_last, '\n');
  if (trailing_lines < rec_lines)
    chunks.back().second = data + prev_start;
  output_cursor = std::distance(data, chunks.back().second);

  return chunks;
}
