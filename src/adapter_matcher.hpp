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

#ifndef SRC_ADAPTER_MATCHER_HPP_
#define SRC_ADAPTER_MATCHER_HPP_

#include "falco_utils.hpp"

#include <array>
#include <cstdint>
#include <string>
#include <thread>  // IWYU pragma: keep
#include <utility>
#include <vector>

struct adapter_matcher {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.05, "pass"},
    std::pair{0.10, "warn"},
    std::pair{std::numeric_limits<double>::max(), "fail"},
  };
  static constexpr auto n_adapters = 6;
  static constexpr auto adapter_size = 12;
  static constexpr auto adapter_names = std::array{
    // clang-format off
    "Illumina Universal Adapter",
    "Illumina Small RNA 3' Adapter",
    "Illumina Small RNA 5' Adapter",
    "Nextera Transposase Sequence",
    "PolyA",
    "PolyG",
    // clang-format on
  };
  static constexpr auto adapters = std::array{
    // clang-format off
    "AGATCGGAAGAG",  // Illumina Universal Adapter
    "TGGAATTCTCGG",  // Illumina Small RNA 3' Adapter
    "GATCGTCGGACT",  // Illumina Small RNA 5' Adapter
    "CTGTCTCTTATA",  // Nextera Transposase Sequence
    "AAAAAAAAAAAA",  // PolyA
    "GGGGGGGGGGGG",  // PolyG
    // clang-format on
  };
  std::array<std::uint64_t, n_adapters> encoded_adapters{};
  std::vector<std::array<std::uint64_t, n_adapters>> adap_counts;

  adapter_matcher();

  [[nodiscard]] static auto
  match_adapter(const auto &t, const auto m, auto adap) -> std::uint32_t {
    static constexpr auto adap_mask = (1ul << adapter_size * nibble_size) - 1ul;
    std::uint64_t t_enc{};
    auto i = 0u;
    for (; i + 1 < adapter_size; ++i)
      t_enc = (t_enc << nibble_size) + encode_nibble(t[i]);
    while (i < m) {
      t_enc = (t_enc << nibble_size) + encode_nibble(t[i++]);
      if (((t_enc ^ adap) & adap_mask) == 0)
        return i - adapter_size;
    }
    return m;
  }

  auto
  match_adapters(const auto seq, const auto len) {
    if (len < adapter_size) [[unlikely]]
      return;
    for (auto i = 0; i < n_adapters; ++i)
      // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-constant-array-index)
      if (const auto p = match_adapter(seq, len, encoded_adapters[i]); p < len)
        ++adap_counts[p][i];
  }

  auto
  resize(const std::uint32_t updated_length) {
    adap_counts.resize(updated_length);
  }

  auto
  operator+=(const adapter_matcher &rhs) -> const adapter_matcher &;

  [[nodiscard]] auto
  string(const std::uint64_t n_reads) const -> std::string;
};

#endif  // SRC_ADAPTER_MATCHER_HPP_
