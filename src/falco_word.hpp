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

#ifndef FALCO_WORD_HPP_
#define FALCO_WORD_HPP_

#include <algorithm>
#include <cstdint>
#include <format>
#include <iterator>
#include <ranges>
#include <string>

struct falco_word {
  static constexpr auto shift_for_width_bits = 56u;
  static constexpr auto bits_per_base = 2u;
  static constexpr auto mask = 3u;
  static constexpr auto max_lo_lim = 32u;
  static constexpr auto max_hi_lim = 56u;
  std::uint64_t lo{};
  std::uint64_t hi{};

  [[nodiscard]] auto
  operator<=>(const falco_word &) const = default;

  falco_word(auto b, const std::uint64_t w) {
    static constexpr auto encode = [](const auto c) {
      return (c >> 1) & mask;  // Ns are counted as G so must be subtracted
    };
    static const auto enc_shift = [&](auto &x, auto &c) {
      x = (x << bits_per_base) | encode(*c++);
    };
    const auto lo_lim = b + (w > max_lo_lim ? max_lo_lim : w);
    const auto hi_lim = b + (w < max_hi_lim ? w : max_hi_lim);
    while (b < lo_lim)
      enc_shift(lo, b);
    while (b < hi_lim)
      enc_shift(hi, b);
    hi |= (w << shift_for_width_bits);
  }

  [[nodiscard]] static auto
  string_impl(auto word, auto n_bases) {
    static constexpr auto bases = "ACTG";
    std::string r;
    for (auto i = 0u; i < n_bases; ++i, word >>= bits_per_base)
      r += bases[word & mask];
    return r;
  }

  [[nodiscard]] auto
  string() const {
    const std::uint32_t n_bases = (hi >> shift_for_width_bits);
    auto r_lo = string_impl(lo, std::min(n_bases, max_lo_lim));
    std::ranges::reverse(r_lo);
    auto r_hi = string_impl(hi, n_bases - std::min(n_bases, max_lo_lim));
    std::ranges::reverse(r_hi);
    return r_lo + r_hi;
  }

  [[nodiscard]] auto
  hash() const noexcept {
    static constexpr auto hashfun = std::hash<std::uint64_t>{};
    const auto h = hashfun(lo);
    return h ^ ((hashfun(hi) + 0x517cc1b727220a95) + (h << 6) + (h >> 2));
  }
};

template <> struct std::formatter<falco_word> : std::formatter<std::string> {
  auto
  format(const falco_word &fw, auto &ctx) const {
    return std::formatter<std::string>::format(fw.string(), ctx);
  }
};

template <> struct std::hash<falco_word> {
  [[nodiscard]] auto
  operator()(const falco_word &fw) const noexcept {
    return fw.hash();
  }
};

#endif  // FALCO_WORD_HPP_
