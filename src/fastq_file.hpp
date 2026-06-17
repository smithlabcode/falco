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

#include "falco_utils.hpp"  // for falco_thread_pool

#ifdef HAVE_ISAL
#include <isa-l/igzip_lib.h>
#endif  // HAVE_ISAL

#include <htslib/bgzf.h>
#include <htslib/sam.h>

#include <fcntl.h>     // for open, O_RDONLY
#include <sys/mman.h>  // for mmap, munmap, MAP_FAILED, MAP_PRIVATE
#include <unistd.h>    // for close, sysconf, _SC_PAGESIZE

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <concepts>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <iterator>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
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
  size() const { return static_cast<std::int32_t>(std::distance(r, o)) - 1; }
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
  auto itr = n + 1;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

  const auto next_newline = [end_itr](auto &itr) {
    itr = std::find(itr, end_itr, '\n');
  };

  // clang-format off
  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto r = itr;

  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto o = itr;

  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto q = itr;

  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto e = itr;
  // clang-format on

  cursor = e;
  return {n, r, o, q, e};
}

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
  if (munmap(static_cast<void *>(buf.data), buf.sz))
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to unmap file");
  buf.data = nullptr;
  buf.sz = 0;
}

struct fastq_file {
  using rec_t = fqrec;
  static constexpr auto min_buf_size = 64 * 1024;
  std::int64_t buf_size{};
  std::int64_t filesize{};
  fastq_buffer buf{};
  std::int64_t start_pos_in_file{};
  std::int64_t stop_pos_in_file{};
  std::int64_t cursor{};
  int fd{};

  fastq_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size < min_buf_size ? min_buf_size : buf_size},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    stop_pos_in_file{buf_size},  // init this way because used as sentinel
    fd{open(std::data(filename), O_RDONLY, 0)} {
    if (fd < 0)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
  }

  // clang-format off
  [[nodiscard]] auto get_cursor() const { return cursor; }
  auto set_cursor(const auto c) { cursor = c; }
  // clang-format on

  // clang-format off
  // delete copy and assignment
  fastq_file(const fastq_file &) = delete;
  auto operator=(const fastq_file &) -> fastq_file & = delete;
  auto operator=(fastq_file &&) noexcept -> fastq_file & = delete;
  // default move for emplace
  fastq_file(fastq_file &&) noexcept = default;
  // clang-format on

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
};

struct fastq_bgzf_file {
  using rec_t = fqrec;
  std::int64_t buf_size{};  // size of allocated buffer
  fastq_buffer buf{};
  std::vector<char> buffer;
  std::int64_t sentinel_position{};
  std::int64_t cursor{};  // position in buffer
  std::unique_ptr<BGZF, int (*)(BGZF *)> f;

  fastq_bgzf_file(const std::string &filename, const std::int64_t buf_size,
                  falco_thread_pool &t) :
    buf_size{buf_size}, buffer(buf_size), sentinel_position{buf_size},
    f(bgzf_open(std::data(filename), "r"), &bgzf_close) {
    static constexpr auto bgzf_fmt_code = 2;  // from bgzf.h
    if (!f)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
    if (t.n_threads() > 0 && bgzf_compression(f.get()) == bgzf_fmt_code) {
      // threads can be used
      const auto r = bgzf_thread_pool(f.get(), t.t.pool, t.t.qsize);
      if (r < 0)
        throw std::runtime_error("failed to set thread pool");
    }
    buf.data = std::data(buffer);
  }

  // clang-format off
  [[nodiscard]] auto get_cursor() const { return cursor; }
  auto set_cursor(const auto c) { cursor = c; }
  // clang-format on

  // clang-format off
  // delete copy and assignment
  fastq_bgzf_file(const fastq_bgzf_file &) = delete;
  auto operator=(const fastq_bgzf_file &) -> fastq_bgzf_file & = delete;
  auto operator=(fastq_bgzf_file &&) noexcept -> fastq_bgzf_file & = delete;
  // default move for emplace
  fastq_bgzf_file(fastq_bgzf_file &&) noexcept = default;
  ~fastq_bgzf_file() = default;
  // clang-format on

  [[nodiscard]] operator bool() const { return sentinel_position == buf_size; }

  auto
  load_next() {
    if (cursor > 0) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      std::copy_n(std::data(buffer) + cursor, buf.sz - cursor,
                  std::data(buffer));
      cursor = buf.sz - cursor;  // backup to after previous data
    }
    const auto n_bytes = buf_size - cursor;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const auto r = bgzf_read(f.get(), std::data(buffer) + cursor, n_bytes);
    if (r < 0) {
      // ADS: cleanup
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed reading gz input file");
    }
    buf.sz = cursor + r;
    sentinel_position = buf.sz;
    cursor = 0;  // cursor always moves to zero because buffer is not mmapped
  }
};

#ifdef HAVE_ISAL

class fastq_gz_file {
  static constexpr auto inflate_err_msg =  //
    R"(Failure during decompression by ISAL. Error code is {}.  Please check that
the input file is not corrupted by decompressing with gunzip. If the input file is
not corrupted please file a bug report at the Falco repo on GitHub.)";

public:
  using rec_t = fqrec;
  fastq_buffer buf;
  std::int64_t cursor{};

private:
  inflate_state state{};
  isal_gzip_header gz_hdr{};
  std::vector<char> outbuf;
  std::vector<char> inbuf;
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> in;

public:
  // clang-format off
  [[nodiscard]] auto get_cursor() const { return cursor; }
  [[nodiscard]] auto get_buf_end() const { return buf.data + buf.sz; }
  auto set_cursor(const auto c) { cursor = c; }
  // clang-format on

  explicit fastq_gz_file(const std::string &filename,
                         const std::int64_t bufsize) :
    outbuf(bufsize / 2), inbuf(bufsize / 2),
    in(std::fopen(std::data(filename), "r"), &std::fclose) {
    if (in == nullptr)
      throw std::runtime_error("failed to open " + filename);

    isal_inflate_init(&state);
    state.crc_flag = ISAL_GZIP_NO_HDR_VER;
    update_state_in();

    // Actually read and save the header info
    isal_gzip_header_init(&gz_hdr);
    const auto ret = isal_read_gzip_header(&state, &gz_hdr);
    if (ret != ISAL_DECOMP_OK)
      throw std::runtime_error("failed to read gz header from: " + filename);
    buf.sz = 0;
    buf.data = std::data(outbuf);
  }

  [[nodiscard]] operator bool() const {
    // ADS: check both conditions below: small files might hit eof on first read
    return !std::feof(in.get()) || state.avail_in > 0;
  }

  auto
  load_next() {
    if (cursor > 0)
      shift_buffer();
    // ADS: below likely won't happen just after reading header in constructor
    if (state.avail_in == 0)
      update_state_in();
    if (state.block_state == ISAL_BLOCK_FINISH && !process_header()) {
      // ADS: we only reach this point if we have multiple gz in the same file
      return;
    }

    // ADS: buf.sz is right after previous decompressed but not-yet-used data
    const std::int64_t remaining_capacity = std::ssize(outbuf) - buf.sz;

    state.avail_out = remaining_capacity;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    state.next_out = reinterpret_cast<std::uint8_t *>(get_buf_end());

    const auto r = isal_inflate(&state);
    if (r != ISAL_DECOMP_OK && r != ISAL_END_INPUT)
      throw std::runtime_error(std::format(inflate_err_msg, r));

    assert(remaining_capacity >= state.avail_out);
    const std::int64_t n_bytes = remaining_capacity - state.avail_out;

    buf.sz += n_bytes;
    cursor = 0;
  }

private:
  auto
  update_state_in() -> void {
    const auto fread = [&](auto &b, const auto n) {
      return static_cast<std::uint32_t>(std::fread(b, 1, n, in.get()));
    };
    // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
    state.next_in = reinterpret_cast<std::uint8_t *>(std::data(inbuf));
    state.avail_in = fread(state.next_in, std::size(inbuf));
  }

  [[nodiscard]] auto
  process_header() -> bool {
    isal_inflate_reset(&state);
    state.crc_flag = ISAL_GZIP_NO_HDR;
    isal_gzip_header_init(&gz_hdr);  // process extra headers
    // skip the next header
    const auto ret = isal_read_gzip_header(&state, &gz_hdr);
    if (ret != ISAL_DECOMP_OK)  // allow trailing junk
      return false;
    return true;
  }

  auto
  shift_buffer() -> void {
    const auto buf_data = std::data(outbuf);
    const auto n_bytes_to_keep = buf.sz - cursor;
    // ADS: need to check the conditions when these might overlap
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    std::memcpy(buf_data, buf_data + cursor, n_bytes_to_keep);
    buf.sz = n_bytes_to_keep;
    cursor = 0;
  }
};

#else  // use bgzf for ordinary gz files

struct fastq_gz_file {
  using rec_t = fqrec;
  std::int64_t buf_size{};  // size of allocated buffer
  std::int64_t buf_used{};
  fastq_buffer buf{};
  std::vector<char> buffer;
  std::int64_t cursor{};  // position in buffer
  std::unique_ptr<BGZF, int (*)(BGZF *)> f;

  fastq_gz_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size}, buf_used{buf_size}, buffer(buf_size),
    f(bgzf_open(std::data(filename), "r"), &bgzf_close) {
    if (!f)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
    buf.data = std::data(buffer);
  }

  // clang-format off
  [[nodiscard]] auto get_cursor() const { return cursor; }
  auto set_cursor(const auto c) { cursor = c; }
  // clang-format on

  // clang-format off
  // delete copy and assignment
  fastq_gz_file(const fastq_gz_file &) = delete;
  auto operator=(const fastq_gz_file &) -> fastq_gz_file & = delete;
  auto operator=(fastq_gz_file &&) noexcept -> fastq_gz_file & = delete;
  // default move for emplace
  fastq_gz_file(fastq_gz_file &&) noexcept = default;
  ~fastq_gz_file() = default;
  // clang-format on

  [[nodiscard]] operator bool() const { return buf_used == buf_size; }

  auto
  load_next() {
    if (cursor > 0) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      std::copy_n(buf.data + cursor, buf.sz - cursor, buf.data);
      cursor = buf.sz - cursor;  // rewind to after previous data
    }
    const auto n_bytes = buf_size - cursor;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const auto r = bgzf_read(f.get(), buf.data + cursor, n_bytes);
    if (r < 0) {
      buf_used = 0;
      return;
    }
    buf.sz = cursor + r;
    buf_used = buf.sz;
    cursor = 0;  // cursor always moves to zero because buffer is not mmapped
  }
};

#endif  // HAVE_ISAL

[[nodiscard]] static inline auto
get_chunks_fastq_impl(auto &fq, const std::int64_t n_chunks) {
  static constexpr auto rec_lines = 4;  // FASTQ
  // clang-format off
  const auto not_read_start = [](const auto s, const auto p) {
    // ADS: could get confused if '+' lines have full name info
    return s[p] != '@' || s[p-1] != '\n' || (s[p-2] == '+' && s[p-3] == '\n');
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
  std::int64_t start_pos = fq.get_cursor();
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
  fq.set_cursor(chunks.back().second);
  return chunks;
}

// specialization to these two classes to distinguish from BAM/SAM
template <typename T>
concept fastq_like =
  std::same_as<T, fastq_file> || std::same_as<T, fastq_bgzf_file> ||
  std::same_as<T, fastq_gz_file>;

template <fastq_like T>
[[nodiscard]] static inline auto
get_chunks(T &fq, const std::int64_t n_chunks) {
  assert(n_chunks > 0);
  const auto add_offset = [d = fq.buf.data](const auto &x) {
    return std::pair{d + x.first, d + x.second};
  };
  return get_chunks_fastq_impl(fq, n_chunks) |
         std::views::transform(add_offset) | std::ranges::to<std::vector>();
}

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

[[nodiscard]] auto
estimate_n_reads_fastq_bgzf(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

[[nodiscard]] auto
estimate_n_reads_fastq_gz(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

#endif  // SRC_FASTQ_FILE_HPP_
