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

#ifndef FALCO_UTILS_HPP_
#define FALCO_UTILS_HPP_

#include <algorithm>
#include <iterator>
#include <numeric>
#include <ranges>

// N (78)10 = (1001110)2
// A (65)10 = (1000001)2
// C (67)10 = (1000011)2
// G (71)10 = (1000111)2
// T (84)10 = (1010100)2

[[nodiscard]] constexpr auto
encode(const char c) {
  return (c >> 1) & 3;  // Ns are counted as G so must be subtracted
}

static constexpr auto nibble_size = 4;

[[nodiscard]] constexpr auto
encode_nibble(const char c) {
  return (c >> 1) & 15;  // N gets separate encoding
}

constexpr auto pct = [](const double a) { return 100.0 * a; };

const auto add = [](auto &a1, auto &a2) {
  std::ranges::transform(a1, a2, std::begin(a1), std::plus{});
};

const auto vec_add = [](auto &v1, const auto &v2) {
  v1.resize(std::max(std::size(v1), std::size(v2)));
  add(v1, v2);
};

const auto two_dim_add = [](auto &v1, const auto &v2) {
  v1.resize(std::max(std::size(v1), std::size(v2)));
  for (auto [a1, a2] : std::views::zip(v1, v2))
    add(a1, a2);
};

[[nodiscard]] inline auto
operator+(auto lhs, const auto &rhs) {
  return lhs += rhs;
}

[[nodiscard]] constexpr auto
as_frac(const auto a, const auto b) {
  return static_cast<double>(a) / static_cast<double>(b);
}

[[nodiscard]] constexpr auto
ipow(const auto b, const auto e) -> std::remove_cvref_t<decltype(b)> {
  return e == 0 ? 1 : (e & 1 ? b : 1) * ipow(b * b, e >> 1);
}

static inline auto
count_nucs(auto seq_itr, const auto sz, auto &tab) {
  const auto seq_end = seq_itr + sz;
  auto out_itr = std::begin(tab);
  while (seq_itr != seq_end)
    ++(*out_itr++)[encode(*seq_itr++)];
}

static inline auto
count_ns(auto seq_itr, const auto sz, auto &tab) {
  const auto seq_end = seq_itr + sz;
  auto out_itr = std::begin(tab);
  while (seq_itr != seq_end)
    *out_itr++ += (*seq_itr++ == 'N');
}

[[nodiscard]] static inline auto
count_gc(auto seq_itr, const auto sz) {
  const auto seq_end = seq_itr + sz;
  auto gc = 0;
  while (seq_itr != seq_end) {
    gc += (*seq_itr != 'N' && (encode(*seq_itr) & 1));
    ++seq_itr;
  }
  return gc;
}

[[nodiscard]] static inline auto
count_quals(auto qual_itr, const auto sz, auto &tab) {
  const auto qual_end = qual_itr + sz;
  auto out_itr = std::begin(tab);
  auto qual_tot = 0;
  while (qual_itr != qual_end) {
    const auto q = *qual_itr++;
    ++(*out_itr++)[q];
    qual_tot += q;
  }
  return qual_tot;
}

static constexpr auto end_module_tag = ">>END_MODULE\n";

[[nodiscard]] static inline auto
tabular_dot(const auto &a) {
  auto total = static_cast<std::remove_cvref_t<decltype(a)>::value_type>(0);
  for (const auto [i, x] : std::views::enumerate(a))
    total += i * x;
  return total;
}

[[nodiscard]] static inline auto
mean_tabular(const auto &a) {
  const auto num = tabular_dot(a);
  const auto denom = std::reduce(std::cbegin(a), std::cend(a));
  return static_cast<double>(num) / static_cast<double>(denom);
}

#endif  // FALCO_UTILS_HPP_
