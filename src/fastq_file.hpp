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

#include <htslib/bgzf.h>
#include <htslib/thread_pool.h>

#include <fcntl.h>
#include <sys/mman.h>
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

// clang-format off
struct fqrec {
  using pos_t = char *;
  pos_t n{};  // start of "name"
  pos_t r{};  // start of "read"
  pos_t o{};  // start of "other"
  pos_t q{};  // start of "quality" scores
  pos_t e{};  // end of the record
  [[nodiscard]] auto
  size() const -> std::uint32_t { return std::distance(r, o) - 1; }
  [[nodiscard]] operator bool() const { return n != nullptr; }
  [[nodiscard]] auto
  string() const -> std::string { return {n, e}; }
};
[[nodiscard]] constexpr auto get_name(const fqrec &rec) { return rec.n; }
[[nodiscard]] constexpr auto get_name_end(const fqrec &rec) { return rec.r - 1; }
[[nodiscard]] constexpr auto get_seq(const fqrec &rec) { return rec.r; }
[[nodiscard]] constexpr auto get_seq_end(const fqrec &rec) { return rec.o - 1; }
[[nodiscard]] constexpr auto get_seq_size(const fqrec &rec) { return std::size(rec); }
[[nodiscard]] constexpr auto get_qual(const fqrec &rec) { return rec.q; }
[[nodiscard]] constexpr auto get_qual_end(const fqrec &rec) { return rec.e - 1; }
[[nodiscard]] constexpr auto get_qual_size(const fqrec &rec) { return std::size(rec); }
// clang-format on

[[nodiscard]] inline auto
get_next(fqrec::pos_t &cursor, const fqrec::pos_t end_itr) -> fqrec {
  // ADS: need to make sure cursor < end_itr or we will move past
  const auto n = cursor;
  auto itr = n + 1;
  itr = std::find(itr, end_itr, '\n');
  if (itr++ == end_itr)
    return {};
  const auto r = itr;

  itr = std::find(itr, end_itr, '\n');
  if (itr++ == end_itr)
    return {};
  const auto o = itr;

  itr = std::find(itr, end_itr, '\n');
  if (itr++ == end_itr)
    return {};
  const auto q = itr;

  itr = std::find(itr, end_itr, '\n');
  if (itr++ == end_itr)
    return {};
  const auto e = itr;

  cursor = e;
  return {n, r, o, q, e};
}

struct fastq_buffer {
  char *data{};       // not necessarily owned
  std::int64_t sz{};  // slight redundancy with vars containing classes
};

[[nodiscard]] inline auto
mmap_fastq(const int fd, const std::int64_t start_pos_in_file,
           const std::int64_t stop_pos_in_file) {
  const auto n_bytes = stop_pos_in_file - start_pos_in_file;
  char *data = static_cast<char *>(
    mmap(nullptr, n_bytes, PROT_READ, MAP_PRIVATE, fd, start_pos_in_file));
  if (data == MAP_FAILED)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to mmap file");
  return fastq_buffer{data, n_bytes};
}

static inline auto
cleanup_mmap_fastq(fastq_buffer &buf) {
  if (buf.data == nullptr)
    return;
  [[maybe_unused]] const int rc = munmap(static_cast<void *>(buf.data), buf.sz);
  buf = {nullptr, 0};
}

struct fastq_file {
  using rec_t = fqrec;
  static constexpr auto min_buf_size = 16 * 4096;
  std::int64_t buf_size{};
  std::int64_t filesize{};
  fastq_buffer buf{};
  std::int64_t start_pos_in_file{};
  std::int64_t stop_pos_in_file{};
  std::int64_t cursor{};
  int fd{};

  fastq_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    stop_pos_in_file{buf_size},  // init this way because used as sentinel
    fd{open(std::data(filename), O_RDONLY, 0)} {
    if (fd < 0)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
    if (buf_size < min_buf_size)
      throw std::runtime_error(
        std::format("Requested buffer size {} smaller than required {}",
                    buf_size, min_buf_size));
  }

  ~fastq_file() {
    if (buf.sz > 0)
      cleanup_mmap_fastq(buf);
    close(fd);  // will always have been opened using a filename
  }

  [[nodiscard]] operator bool() const {
    return stop_pos_in_file - start_pos_in_file == buf_size;
  }

  auto
  load_next() {
    // memory mapped data is page aligned; the data we need is not
    static const auto page_mask = sysconf(_SC_PAGESIZE) - 1;
    std::tie(start_pos_in_file, cursor) = [&] {
      const auto cursor_pos_in_file = start_pos_in_file + cursor;
      return std::tuple(cursor_pos_in_file & (~page_mask),
                        cursor_pos_in_file & page_mask);
    }();
    stop_pos_in_file = std::min(filesize, start_pos_in_file + buf_size);
    if (buf.sz > 0)
      cleanup_mmap_fastq(buf);
    buf = mmap_fastq(fd, start_pos_in_file, stop_pos_in_file);
  }
};

struct fastq_gz_file {
  using rec_t = fqrec;
  std::int64_t buf_size{};  // size of allocated buffer
  std::int64_t filesize{};
  fastq_buffer buf{};
  std::int64_t start_pos_in_file{};  // file offset for buffer start
  std::int64_t stop_pos_in_file{};   // file offset for buffer stop
  std::int64_t cursor{};             // position in buffer
  falco_thread_pool t;
  std::unique_ptr<BGZF, int (*)(BGZF *)> f;

  fastq_gz_file(const std::string &filename, const std::int64_t buf_size,
                const std::uint32_t n_threads = 1) :
    buf_size{buf_size},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    stop_pos_in_file{buf_size}, t(n_threads),
    f(bgzf_open(std::data(filename), "r"), &bgzf_close) {
    static constexpr auto bgzf_fmt_code = 2;  // from bgzf.h
    if (!f)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
    if (n_threads > 1 && bgzf_compression(f.get()) == bgzf_fmt_code) {
      // threads can be used
      const auto r = bgzf_thread_pool(f.get(), t.t.pool, t.t.qsize);
      if (r < 0)
        throw std::runtime_error("failed to set thread pool");
    }
    buf.data = new char[buf_size];
  }

  ~fastq_gz_file() { delete[] buf.data; }

  [[nodiscard]] operator bool() const { return stop_pos_in_file == buf_size; }

  auto
  load_next() {
    if (cursor > 0) {
      std::memmove(buf.data, buf.data + cursor, buf.sz - cursor);
      cursor = buf.sz - cursor;  // backup to after previous data
    }
    start_pos_in_file = cursor;
    const auto n_bytes = buf_size - start_pos_in_file;
    const auto r = bgzf_read(f.get(), buf.data + start_pos_in_file, n_bytes);
    if (r < 0) {
      // ADS: cleanup
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed reading gz input file");
    }
    stop_pos_in_file = start_pos_in_file + r;
    buf.sz = stop_pos_in_file;
    cursor = 0;  // cursor always moves to zero because buffer is not mmapped
  }
};

[[nodiscard]] static inline auto
get_chunks_fastq_impl(auto &fq, const std::int64_t n_chunks)
  -> std::vector<std::pair<std::int64_t, std::int64_t>> {
  static constexpr auto rec_lines = 4;  // FASTQ
  // clang-format off
  const auto not_read_start = [](const auto s, const auto p) {
    return s[p] != '@' || s[p-1] != '\n' || (s[p-2] == '+' && s[p-3] == '\n');
  };
  const auto fwd_to_read_start = [&](const fastq_buffer &buf, auto x) {
    if (x == 0) return x;
    while (x < buf.sz && not_read_start(buf.data, x)) ++x;
    return x;
  };
  const auto rev_to_read_start = [&](const fastq_buffer &buf, auto x) {
    while (x > 0 && (x == buf.sz || not_read_start(buf.data, x))) --x;
    return x;
  };
  // clang-format on
  const auto buf = fq.buf;
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
  if (std::count(buf.data + prev_start, buf.data + buf.sz, '\n') < rec_lines)
    chunks.back().second = prev_start;
  fq.cursor = chunks.back().second;
  return chunks;
}

// specialization to these two classes
template <class T>
  requires std::same_as<T, fastq_file> || std::same_as<T, fastq_gz_file>
[[nodiscard]] static inline auto
get_chunks(T &fq, const std::int64_t n_chunks)
  -> std::vector<std::pair<fqrec::pos_t, fqrec::pos_t>> {
  // const auto chunks = get_chunks_fastq_impl(fq, n_chunks);
  return get_chunks_fastq_impl(fq, n_chunks) |
         std::views::transform([&](const auto &x) {
           return std::pair{fq.buf.data + x.first, fq.buf.data + x.second};
         }) |
         std::ranges::to<std::vector>();
}

#endif  // SRC_FASTQ_FILE_HPP_
