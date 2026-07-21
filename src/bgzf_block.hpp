// SPDX-License-Identifier: MIT; (c) 2026 Andrew D Smith

#ifndef SRC_BGZF_BLOCK_HPP_
#define SRC_BGZF_BLOCK_HPP_

#include <libdeflate.h>

#include <array>
#include <cstdint>
#include <vector>

static constexpr auto max_bgzf_block_size = 65536;

struct bgzf_block_t {
  std::int32_t size{};
  char *out_itr{};
  char *in_itr;

  bgzf_block_t(const std::int32_t size, char *out_itr, char *in_itr) noexcept :
    size{size}, out_itr{out_itr}, in_itr{in_itr} {}
  // clang-format off
  bgzf_block_t(const bgzf_block_t &src) = delete;
  auto operator=(const bgzf_block_t &src) -> bgzf_block_t & = delete;

  bgzf_block_t() = default;
  bgzf_block_t(bgzf_block_t &&src) noexcept = default;
  auto operator=(bgzf_block_t &&src) noexcept -> bgzf_block_t & = default;

  [[nodiscard]] auto data() const -> const char * { return in_itr; }
  [[nodiscard]] operator bool() const { return size > 0; }
  [[nodiscard]] auto operator<=>(const bgzf_block_t &other) const = default;
  // clang-format on

  auto
  decompress() -> void;
};

auto
decompress(auto &&x) -> void {  // free function
  x.decompress();
}

#endif  // SRC_BGZF_BLOCK_HPP_
