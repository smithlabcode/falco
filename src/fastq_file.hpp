// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FASTQ_FILE_HPP_
#define SRC_FASTQ_FILE_HPP_

#include "fqrec.hpp"
#include "task_queue.hpp"

#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iterator>
#include <ranges>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

struct fastq_buffer {
  char *data{};       // not necessarily owned
  std::int64_t sz{};  // slight redundancy with vars containing classes
};

static inline auto
mmap_fastq(const int fd, const std::int64_t start_pos_in_file,
           const std::int64_t stop_pos_in_file, fastq_buffer &buf) {
  const auto n_bytes = stop_pos_in_file - start_pos_in_file;
  char *data = static_cast<char *>(
    mmap(nullptr, n_bytes, PROT_READ, MAP_PRIVATE, fd, start_pos_in_file));
  if (data == MAP_FAILED)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to mmap file");
  buf.data = data;
  buf.sz = n_bytes;
}

static inline auto
cleanup_mmap_fastq(fastq_buffer &buf) {
  if (buf.data == nullptr)
    return;
  [[maybe_unused]] const auto r = munmap(static_cast<void *>(buf.data), buf.sz);
  buf.data = nullptr;
  buf.sz = 0;
}

struct fastq_file {
  static constexpr auto min_buf_size = 64 * 1024;
  std::int64_t buf_size{};
  std::int64_t filesize{};
  fastq_buffer buf{};
  std::int64_t start_pos_in_file{};
  std::int64_t stop_pos_in_file{};
  std::int64_t cursor{};
  int fd{};

  fastq_file(const std::string &filename, const std::int64_t buf_size_arg) :
    buf_size{buf_size_arg + min_buf_size},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    stop_pos_in_file{buf_size_arg + min_buf_size},  // init to use as sentinel
    fd{open(std::data(filename), O_RDONLY, 0)} {
    if (fd < 0)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
  }

  // clang-format off
  // delete copy and assignment
  fastq_file(const fastq_file &) = delete;
  auto operator=(const fastq_file &) -> fastq_file & = delete;
  auto operator=(fastq_file &&) noexcept -> fastq_file & = delete;
  // clang-format on

  fastq_file(fastq_file &&src) noexcept :
    buf_size{src.buf_size},                    //
    filesize{src.filesize},                    //
    buf{src.buf},                              //
    start_pos_in_file{src.start_pos_in_file},  //
    stop_pos_in_file{src.stop_pos_in_file},    //
    cursor{src.cursor},                        //
    fd{dup(src.fd)}                            // <- LOOK
  {}

  auto
  reset() -> void {
    if (buf.sz > 0)
      cleanup_mmap_fastq(buf);
  }

  ~fastq_file() {
    reset();
    close(fd);  // will always have been opened using a filename
  }

  [[nodiscard]] operator bool() const {
    return stop_pos_in_file - start_pos_in_file == buf_size;
  }

  auto
  load_next() -> void {
    // memory mapped data is page aligned but the data we need is not
    static const auto page_mask = sysconf(_SC_PAGESIZE) - 1;
    std::tie(start_pos_in_file, cursor) = [&] {
      const auto pos_in_file = start_pos_in_file + cursor;
      return std::tuple(pos_in_file & (~page_mask), pos_in_file & page_mask);
    }();
    stop_pos_in_file = std::min(filesize, start_pos_in_file + buf_size);
    if (buf.sz > 0)
      cleanup_mmap_fastq(buf);
    mmap_fastq(fd, start_pos_in_file, stop_pos_in_file, buf);
  }

  [[nodiscard]] auto
  get_chunks(const std::int64_t n_chunks) -> fq_chunks_t;
};

[[nodiscard]] static inline auto
get_chunks_fastq_impl(auto &fq, const std::int64_t n_chunks) {
  static constexpr auto rec_lines = 4;  // FASTQ
  // clang-format off
  const auto not_read_start = [](const auto s, const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return s[p] != '@' || (p > 0 && s[p-1] != '\n') ||
      (p > 2 && s[p-2] == '+' && s[p-3] == '\n');
  };
  const auto fwd_to_read_start = [&](const auto &buf, auto pos) {
    if (pos == 0) return pos;
    while (pos < buf.sz && not_read_start(buf.data, pos)) ++pos;
    return pos;
  };
  const auto rev_to_read_start = [&](const auto &buf, auto pos) {
    while (pos > 0 && (pos == buf.sz || not_read_start(buf.data, pos))) --pos;
    return pos;
  };
  // clang-format on
  const auto &buf = fq.buf;  // non-copyable
  const auto n_bytes_available = buf.sz - fq.cursor;
  const auto [chunk_size, remainder] = std::div(n_bytes_available, n_chunks);
  std::vector<std::pair<std::int64_t, std::int64_t>> chunks(n_chunks);
  std::int64_t start_pos = fq.cursor;
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto chunk_beg = fwd_to_read_start(buf, start_pos);
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    const auto chunk_end = fwd_to_read_start(buf, stop_pos);
    chunks[chunk_idx] = {chunk_beg, chunk_end};
    start_pos = stop_pos;
  }
  // make sure final chunk includes only full records
  const auto prev_start = rev_to_read_start(buf, chunks.back().second);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  if (std::count(buf.data + prev_start, buf.data + buf.sz, '\n') < rec_lines)
    chunks.back().second = prev_start;
  fq.cursor = chunks.back().second;
  return chunks;
}

[[nodiscard]] inline auto
fastq_file::get_chunks(const std::int64_t n_chunks) -> fq_chunks_t {
  assert(n_chunks > 0);
  const auto add_offset = [d = buf.data](const auto &x) {
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    return std::pair{d + x.first, d + x.second};
  };
  return get_chunks_fastq_impl(*this, n_chunks) |
         std::views::transform(add_offset) | std::ranges::to<std::vector>();
}

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

inline auto
make_tasks(fastq_file &reads_file,       //
           const std::int64_t n_chunks,  //
           const std::int32_t file_id,   //
           task_queue &tq,               //
           std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes tasks
  reads_file.load_next();
  auto chunks = reads_file.get_chunks(n_chunks);
  for (const auto &chunk : chunks) {
    ++n_tasks;  // ADS: increment before submitting
    tq.push(file_id, fq_task_t(chunk.first, chunk.second));
  }
}

inline auto
reset(fastq_file &reads_file) -> void {
  reads_file.reset();
}

#endif  // SRC_FASTQ_FILE_HPP_
