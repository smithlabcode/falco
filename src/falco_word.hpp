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

#ifndef SRC_FALCO_WORD_HPP_
#define SRC_FALCO_WORD_HPP_

#include <algorithm>
#include <cstdint>
#include <format>
#include <iterator>
#include <ranges>
#include <string>

// conversion for kmers to include 'N'
// N (78)10 = (1001110)2 => 100
// A (65)10 = (1000001)2 => 011
// C (67)10 = (1000011)2 => 010
// G (71)10 = (1000111)2 => 000
// T (84)10 = (1010100)2 => 001

struct falco_word {
  static constexpr auto shift_for_width_bits = 56u;
  static constexpr auto width_bits_removal_mask = 0xffffffffffffff;
  static constexpr auto max_lo_lim = 27u;
  static constexpr auto max_hi_lim = 50u;  // ADS: we could go to 51 bp
  std::uint64_t lo{};
  std::uint64_t hi{};

  [[nodiscard]] auto
  operator<=>(const falco_word &) const = default;

  falco_word(auto b, std::uint64_t w) {
    static constexpr auto fw_encode = [](const auto c) {
      // NOLINTNEXTLINE (cppcoreguidelines-avoid-magic-numbers)
      return ((c >> 1) & 7) ^ 3;
    };
    static const auto enc_shift = [&](auto &x, auto &c) {
      // NOLINTNEXTLINE (cppcoreguidelines-avoid-magic-numbers)
      x = (x * 5) + fw_encode(*c++);
    };
    w = w < max_hi_lim ? w : max_hi_lim;
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
    static constexpr auto extended_bases = "GTCAN";
    static constexpr auto extended_alpha_size = 5;
    std::string r;
    for (auto i = 0u; i < n_bases; ++i) {
      // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-pointer-arithmetic)
      r += extended_bases[word % extended_alpha_size];
      word /= extended_alpha_size;
    }
    return r;
  }

  [[nodiscard]] auto
  string() const {
    const std::uint32_t n_bases = (hi >> shift_for_width_bits);
    auto r_lo = string_impl(lo, std::min(n_bases, max_lo_lim));
    std::ranges::reverse(r_lo);
    const auto hi_clean = hi & width_bits_removal_mask;
    auto r_hi = string_impl(hi_clean, n_bases - std::min(n_bases, max_lo_lim));
    std::ranges::reverse(r_hi);
    return r_lo + r_hi;
  }

  [[nodiscard]] auto
  hash() const noexcept {
    // ADS: from boost multiprecision hash, but for 64 bits
    static constexpr auto magic = 0x517cc1b727220a95;
    static constexpr auto hashfun = std::hash<std::uint64_t>{};
    const auto h = hashfun(lo);
    // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-pointer-arithmetic)
    return h ^ ((hashfun(hi) + magic) + (h << 6) + (h >> 2));
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

[[nodiscard]] inline auto
hash_value(const falco_word &fw) noexcept -> std::uint64_t {
  /// Some containers might need this instead of std::hash for ADL
  return fw.hash();
}

#endif  // SRC_FALCO_WORD_HPP_
