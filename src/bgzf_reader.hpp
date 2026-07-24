// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BGZF_READER_HPP_
#define SRC_BGZF_READER_HPP_

#include "bgzf_block.hpp"

#include <cstdint>
#include <cstdio>
#include <iterator>
#include <memory>
#include <string>

static constexpr auto gzip_header_size = 18;

struct gzip_header {
  static constexpr auto magic1 = 0x1F;
  static constexpr auto magic2 = 0x8B;

  std::uint8_t id1{};       // 0
  std::uint8_t id2{};       // 1
  std::uint8_t cm_eight{};  // 2
  std::uint8_t flg{};       // 3
  std::uint32_t mtime{};    // 7 [4]
  std::uint8_t xfl{};       // 8
  std::uint8_t os{};        // 9
  std::uint16_t xlen{};     // 11 [2]
  char b{};                 // 12
  char c{};                 // 13
  std::uint16_t two{};      // 14 [2]
  std::uint16_t size{};     // 16 [2]
  // 18

  [[nodiscard]] auto
  check_magic() const -> bool {
    return id1 == magic1 && id2 == magic2;  // ADS: check b and c also
  }
};

// reads data and provides serialized compressed chunks for deflation
class bgzf_reader {
private:
  static constexpr auto inbuf_size = 16 * max_bgzf_block_size;

  std::unique_ptr<std::FILE, int (*)(std::FILE *)> fp;
  std::uint64_t filesize{};
  gzip_header gh;
  std::unique_ptr<char[]> inbuf;
  std::unique_ptr<char[]> outbuf;
  char *next_in{};
  char *end_in{};
  char *next_out{};
  char *end_out{};

public:
  bgzf_reader(const std::string &filename, const std::int32_t buf_size);

  operator bool() const { return can_produce_data(); }

  [[nodiscard]] auto
  get_decomp_task(char *out_itr) -> bgzf_block_t;

  [[nodiscard]] auto
  read_data() -> bool;

  auto
  reset() {
    next_out = outbuf.get();
  }

  [[nodiscard]] auto
  task_ready() const -> bool {
    return can_produce_data() && has_out();
  }

private:
  [[nodiscard]] auto
  at_eof() const -> bool {
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
