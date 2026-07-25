// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FASTQ_GZ_FILE_HPP_
#define SRC_FASTQ_GZ_FILE_HPP_

#include "fqrec.hpp"
#include "task_queue.hpp"

#ifdef HAVE_ISAL
#include <isa-l/igzip_lib.h>
#else
#include <htslib/bgzf.h>
#endif  // HAVE_ISAL

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <format>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#ifdef HAVE_ISAL

class fastq_gz_file {
  static constexpr auto min_buf_size = 64 * 1024;
  static constexpr auto inflate_err_msg =  //
    R"(Failure during decompression by ISAL. Error code is {}.  Please check that the
input file is not corrupted by decompressing with gunzip. If the input file is
not corrupted please file a bug report at the Falco repo on GitHub.)";

public:
  using rec_t = fqrec;

private:
  char *buf_data{};       // start of the outbuf
  std::int64_t buf_sz{};  // amount of output buffer used
  std::int64_t cursor{};
  inflate_state state{};
  isal_gzip_header gz_hdr{};
  std::vector<char> outbuf;
  std::vector<char> inbuf;
  bool isal_ok{true};
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> in;

public:
  // clang-format off
  [[nodiscard]] auto get_cursor() const -> std::int64_t { return cursor; }
  auto set_cursor(const auto c) { cursor = c; }
  // clang-format on
  [[nodiscard]] auto
  get_buf_end() const {
    // NOLINTNEXTLINE(*-reinterpret-cast,*-pointer-arithmetic)
    return reinterpret_cast<std::uint8_t *>(buf_data + buf_sz);
  }

  explicit fastq_gz_file(const std::string &filename,
                         const std::int64_t buf_size) :
    outbuf(buf_size / 2), inbuf(buf_size / 2),
    in(std::fopen(std::data(filename), "r"), &std::fclose) {
    if (in == nullptr)
      throw std::runtime_error("failed to open " + filename);
    if (buf_size < min_buf_size)
      throw std::runtime_error(std::format(
        "requested buffer too small {} (min is {})", buf_size, min_buf_size));

    isal_gzip_header_init(&gz_hdr);         // header
    isal_inflate_init(&state);              // inflate
    state.crc_flag = ISAL_GZIP_NO_HDR_VER;  // CRC

    read_data();

    const auto ret = isal_read_gzip_header(&state, &gz_hdr);
    if (ret != ISAL_DECOMP_OK)
      throw std::runtime_error("failed to read gz header from: " + filename);

    buf_sz = 0;
    buf_data = std::data(outbuf);
  }

  [[nodiscard]] operator bool() const {  // does compressed data remain?
    return isal_ok && (!std::feof(in.get()) || state.avail_in > 0);
  }

  auto
  load_next() -> void {
    if (cursor > 0)
      shift_output_buffer();
    while (true) {
      if (need_data())
        read_data();
      if (state.block_state == ISAL_BLOCK_FINISH && !next_gzip_header_ok()) {
        isal_ok = false;  // end of file, but next data not valid file start
        break;
      }
      const auto prev_avail_out = update_output_state();
      if (prev_avail_out == 0)
        break;
      if (const auto r = isal_inflate(&state); r != ISAL_DECOMP_OK)
        throw std::runtime_error(std::format(inflate_err_msg, r));
      buf_sz += (prev_avail_out - state.avail_out);
      if (work_done())
        break;
    }
  }

  [[nodiscard]] auto
  get_chunks(const std::int64_t n_chunks) -> fq_chunks_t;

private:
  [[nodiscard]] auto
  check_gzip_magic() const -> bool {
    static constexpr auto magic1 = 31;
    static constexpr auto magic2 = 139;
    return state.avail_in > 1 && state.next_in[0] == magic1 &&
           state.next_in[1] == magic2;
  }

  [[nodiscard]] auto
  update_output_state() -> std::int64_t {
    const std::int64_t available = std::ssize(outbuf) - buf_sz;
    state.avail_out = available;
    state.next_out = get_buf_end();
    return available;
  }

  [[nodiscard]] auto
  work_done() const -> bool {
    return state.avail_in == 0 && std::feof(in.get());
  }

  [[nodiscard]] auto
  need_data() const -> bool {
    return state.avail_in == 0 && !std::feof(in.get());
  }

  auto
  read_data() -> void {
    const auto fread =
      [&](auto &b, const auto n) {  // cppcheck-suppress constParameterReference
        return static_cast<std::uint32_t>(std::fread(b, 1, n, in.get()));
      };
    // NOLINTNEXTLINE(*-reinterpret-cast)
    state.next_in = reinterpret_cast<std::uint8_t *>(std::data(inbuf));
    state.avail_in = fread(state.next_in, std::size(inbuf));
  }

  [[nodiscard]] auto
  next_gzip_header_ok() -> bool {
    if (check_gzip_magic()) {
      isal_inflate_reset(&state);
      state.crc_flag = ISAL_GZIP;  // process extra headers
      return true;
    }
    return false;
  }

  auto
  shift_output_buffer() -> void {
    // const auto buf_data = std::data(outbuf);
    const auto n_bytes_to_keep = buf_sz - cursor;
    assert(n_bytes_to_keep < cursor);  // ADS: memcpy breaks if overlap
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    std::memcpy(buf_data, buf_data + cursor, n_bytes_to_keep);
    buf_sz = n_bytes_to_keep;
    cursor = 0;
  }
};

#else  // use bgzf for ordinary gz files

struct fastq_gz_file {
  using rec_t = fqrec;
  static constexpr auto min_buf_size = 64 * 1024;
  char *buf_data{};         // start of the outbuf
  std::int64_t buf_sz{};    // amount of output buffer used
  std::int64_t buf_size{};  // size of allocated buffer
  std::int64_t buf_used{};
  std::vector<char> outbuf;
  std::int64_t cursor{};  // position in buffer
  std::unique_ptr<BGZF, int (*)(BGZF *)> f;

  fastq_gz_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size}, buf_used{buf_size}, outbuf(buf_size),
    f(bgzf_open(std::data(filename), "r"), &bgzf_close) {
    if (!f)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
    if (buf_size < min_buf_size)
      throw std::runtime_error(std::format(
        "requested buffer too small {} (min is {})", buf_size, min_buf_size));
    buf_data = std::data(outbuf);
  }

  // clang-format off
  [[nodiscard]] auto get_cursor() const -> std::int64_t { return cursor; }
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
  load_next() -> void {
    if (cursor > 0) {
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
      std::copy_n(buf_data + cursor, buf_sz - cursor, buf_data);
      cursor = buf_sz - cursor;  // rewind to after previous data
    }
    const auto n_bytes = buf_size - cursor;
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    const auto r = bgzf_read(f.get(), buf_data + cursor, n_bytes);
    if (r < 0) {
      buf_used = 0;
      return;
    }
    buf_sz = cursor + r;
    buf_used = buf_sz;
    cursor = 0;  // cursor always moves to zero
  }

  [[nodiscard]] auto
  get_chunks(const std::int64_t n_chunks) -> fq_chunks_t;
};

#endif  // HAVE_ISAL

[[nodiscard]] auto
estimate_n_reads_fastq_gz(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

inline auto
make_tasks(fastq_gz_file &reads_file,    //
           const std::int64_t n_chunks,  //
           const std::int32_t file_id,   //
           task_queue &tq,               //
           std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes tasks
  reads_file.load_next();
  auto chunks = reads_file.get_chunks(n_chunks);
  for (const auto &chunk : chunks) {
    ++n_tasks;
    tq.push(file_id, fq_task_t(chunk.first, chunk.second));
  }
}

#endif  // SRC_FASTQ_GZ_FILE_HPP_
