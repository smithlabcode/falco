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

#include "adapter_matcher.hpp"
#include "falco_utils.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <format>
#include <ranges>
#include <string>
#include <vector>

adapter_matcher::adapter_matcher(const std::uint32_t max_read_len) :
  max_read_len{max_read_len},
  adapter_counts(max_read_len, std::array<std::uint64_t, n_adapters>{}) {
  const auto encode_adap = [&](const auto &a) {
    std::uint64_t x{};
    for (const auto i : std::views::iota(0, adapter_size))
      x = (x << nibble_size) + encode_nibble(a[i]);
    return x;
  };
  for (const auto i : std::views::iota(0, n_adapters))
    encoded_adapters[i] = encode_adap(adapters[i]);
}

auto
adapter_matcher::operator+=(const adapter_matcher &rhs)
  -> const adapter_matcher & {
  two_dim_add(adapter_counts, rhs.adapter_counts);
  max_read_len = std::max(max_read_len, rhs.max_read_len);
  return *this;
}

[[nodiscard]] auto
adapter_matcher::string(const std::uint64_t n_reads) const -> std::string {
  constexpr auto as_frac = [](const auto a, const auto b) {
    return static_cast<double>(a) / static_cast<double>(b);
  };
  auto r = std::format(">>Adapter Content\t{}\n", "pass");
  for (auto i = 0; i < max_read_len; ++i) {
    r += std::format("{}", i + 1);
    for (const auto j : adapter_counts[i])
      r += std::format("\t{:.6g}", as_frac(j, n_reads));
    r += '\n';
  }
  return r + end_module_tag;
}
