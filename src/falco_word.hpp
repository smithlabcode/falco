// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_WORD_HPP_
#define SRC_FALCO_WORD_HPP_

#include <algorithm>
#include <cstdint>
#include <format>
#include <iterator>
#include <ranges>
#include <span>
#include <string>

// conversion for kmers to include 'N'
// N (78)10 = (1001110)2 => 100
// A (65)10 = (1000001)2 => 011
// C (67)10 = (1000011)2 => 010
// G (71)10 = (1000111)2 => 000
// T (84)10 = (1010100)2 => 001

struct falco_word {
  static constexpr std::span extended_bases = "GTCAN";
  static constexpr auto extended_alpha_size = 5;
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
      return ((c >> 1) & 7) ^ 3;  // NOLINT(*-avoid-magic-numbers)
    };
    static const auto enc_shift = [&](auto &x, auto &c) {
      x = (x * extended_alpha_size) + fw_encode(*c++);
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
    std::string r;
    for (auto i = 0u; i < n_bases; ++i) {
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
    // ADS: this might already be in the boost unordered_map header
    static constexpr auto left_shift = 6;
    static constexpr auto magic = 0x517cc1b727220a95;
    static constexpr auto hashfun = std::hash<std::uint64_t>{};
    const auto h = hashfun(lo);
    return h ^ ((hashfun(hi) + magic) + (h << left_shift) + (h >> 2));
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
