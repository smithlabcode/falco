// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BGZF_BLOCK_HPP_
#define SRC_BGZF_BLOCK_HPP_

#include <cstdint>
#include <memory>
#include <vector>

static constexpr auto max_bgzf_block_size = 65536;

struct bgzf_block_t {
  using iterator = char *;

  // ADS: because out_itr and in_itr reference separate blocks, either we need
  // to track one of the sizes or we would have to add an end iterator for one
  // of them.
  std::int32_t size{};
  iterator out_itr{};
  iterator in_itr{};
  bgzf_block_t(const std::int32_t size,
               const iterator out_itr,
               const iterator in_itr) noexcept :
    size{size}, out_itr{out_itr}, in_itr{in_itr} {}

  // clang-format off
  bgzf_block_t(const bgzf_block_t &src) = delete;
  auto operator=(const bgzf_block_t &src) -> bgzf_block_t & = delete;
  bgzf_block_t() = default;
  ~bgzf_block_t() = default;
  bgzf_block_t(bgzf_block_t &&src) noexcept = default;
  auto operator=(bgzf_block_t &&src) noexcept -> bgzf_block_t & = default;
  operator bool() const { return size > 0; }
  // clang-format on

  auto
  decompress() -> void;
};

inline auto
decompress(bgzf_block_t &x) -> void {  // free function
  x.decompress();
}

#endif  // SRC_BGZF_BLOCK_HPP_
