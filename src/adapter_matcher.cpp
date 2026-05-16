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
#include <functional>
#include <iterator>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

adapter_matcher::adapter_matcher() {
  const auto encode_adap = [&](const auto &a) {
    std::uint64_t x{};
    for (const auto i : std::views::iota(0, adapter_size))
      // cppcheck-suppress useStlAlgorithm
      x = (x << nibble_size) + encode_nibble(a[i]);
    return x;
  };
  // NOLINTBEGIN (*-pro-bounds-constant-array-index)
  for (const auto i : std::views::iota(0, n_adapters))
    // cppcheck-suppress useStlAlgorithm
    encoded_adapters[i] = encode_adap(adapters[i]);
  // NOLINTEND (*-pro-bounds-constant-array-index)
}

auto
adapter_matcher::operator+=(const adapter_matcher &rhs)
  -> const adapter_matcher & {
  two_dim_add(adap_counts, rhs.adap_counts);
  return *this;
}

[[nodiscard]] auto
adapter_matcher::string(const std::uint64_t n_reads,
                        const std::uint32_t n_pos) const -> std::string {
  static constexpr auto start_module_tag = ">>Adapter Content\t{}\n";
  static constexpr auto header = "#Position\t";
  static const auto to_flat = [](const auto &fmt, const auto data) {
    return data | std::views::transform(fmt) | std::views::join |
           std::ranges::to<std::string>();
  };

  // need grade first to format module start
  const auto mc =
    std::ranges::max(adap_counts | std::views::transform(std::ranges::max));
  const auto mcf = as_frac(mc, n_reads);
  const auto grade = get_grade(grade_cutoffs, mcf);

  auto r = std::format(start_module_tag, grade);
  r += header + to_flat([](const auto x) { return std::format("\t{}", x); },
                        adapter_names);
  r += '\n';

  auto cumulative = adap_counts;
  for (auto [prev, curr] : cumulative | std::views::pairwise)
    std::ranges::transform(curr, prev, std::begin(curr), std::plus{});

  const auto lim =
    n_pos < std::size(adap_counts) ? n_pos : std::size(adap_counts);
  cumulative.resize(lim);
  const auto fmt_pct_of_reads = [n_reads](const auto c) {
    return std::format("\t{:.6f}", pct(as_frac(c, n_reads)));
  };
  for (const auto [idx, cumul] : std::views::enumerate(cumulative))
    r += std::format("{}{}\n", idx + 1, to_flat(fmt_pct_of_reads, cumul));

  return r + end_module_tag;
}
