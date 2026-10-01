// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BGZF_READER_HPP_
#define SRC_BGZF_READER_HPP_

#include "bgzf_block.hpp"

#include <cstdint>
#include <cstdio>
#include <iterator>
#include <memory>
#include <string>

// reads data and provides serialized compressed chunks for deflation
class bgzf_reader {
private:
  static constexpr auto inbuf_size = 16 * max_bgzf_block_size;

  std::unique_ptr<std::FILE, int (*)(std::FILE *)> fp;
  std::uint64_t filesize{};
  std::unique_ptr<char[]> inbuf;   // NOLINT(cppcoreguidelines-avoid-c-arrays)
  std::unique_ptr<char[]> outbuf;  // NOLINT(cppcoreguidelines-avoid-c-arrays)
  char *next_in{};
  char *end_in{};
  char *next_out{};
  char *end_out{};

public:
  bgzf_reader(const std::string &filename, const std::int64_t buf_size);

  operator bool() const { return can_produce_data(); }

  [[nodiscard]] auto
  get_decomp_task(char *out_itr) -> bgzf_block_t;

  [[nodiscard]] auto
  read_data() -> bool;

  auto
  release() {
    // ADS: intended analogous to monotonic_buffer_resource release
    next_out = outbuf.get();
  }

  [[nodiscard]] auto
  task_ready() const -> bool {
    return can_produce_data() && has_out();
  }

  auto
  reset() {
    inbuf.reset(nullptr);
    outbuf.reset(nullptr);
  }

private:
  [[nodiscard]] auto
  at_eof() const -> bool {
    // ADS: directly checking eof is weird because the block structure of the
    // file means we should read the final byte of the file without hitting eof.
    const auto ft = std::ftell(fp.get());
    return static_cast<std::uint64_t>(ft) == filesize || ft < 0;
  }

  [[nodiscard]] auto
  can_produce_data() const -> bool {
    return !at_eof() || next_in != end_in;
  }

  [[nodiscard]] auto
  has_out() const -> bool {
    return std::distance(next_out, end_out) >= max_bgzf_block_size;
  }
};

#endif  // SRC_BGZF_READER_HPP_
