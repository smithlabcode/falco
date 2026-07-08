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

#ifndef SRC_FALCO_UTILS_HPP_
#define SRC_FALCO_UTILS_HPP_

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <concepts>
#include <cstdint>
#include <functional>
#include <iterator>
#include <numeric>
#include <ranges>
#include <string>
#include <type_traits>
#include <vector>

static constexpr auto cytosine_index = 3;
static constexpr auto guanine_index = 1;

namespace falco {
static constexpr auto alphabet_size = 4;
static constexpr auto gc_content_arra_size = 101;
using nuc_array = std::array<std::uint64_t, alphabet_size>;
using gc_content_array = std::array<std::uint64_t, gc_content_arra_size>;
}  // namespace falco

static constexpr std::int64_t gigabytes = 1024 * 1024 * 1024;
static constexpr std::int64_t megabytes = 1024 * 1024;
static constexpr std::int64_t kilobytes = 1024;

[[nodiscard]] auto
smooth_gc_content(const falco::gc_content_array &data,
                  const std::int64_t window_size) -> std::vector<double>;

[[nodiscard]] auto
get_theoretical_distribution(const falco::gc_content_array &gc,
                             const std::uint64_t total_count)
  -> std::vector<double>;

[[nodiscard]] auto
sum_deviation_from_normal(const falco::gc_content_array &gc) -> double;

[[nodiscard]] inline constexpr auto
duration(const auto start, const auto stop) {
  const auto d = stop - start;
  // ADS: 'count()' because macos has locale issues formatting times
  return std::chrono::duration_cast<std::chrono::duration<double>>(d).count();
};

inline constexpr auto end_module_tag = ">>END_MODULE\n";

// N (78)10 = (1001110)2
// A (65)10 = (1000001)2
// C (67)10 = (1000011)2
// G (71)10 = (1000111)2
// T (84)10 = (1010100)2

[[nodiscard]] inline constexpr auto
encode(const char c) {
  return (c >> 1) & 3;  // Ns are counted as G so must be subtracted
}

[[nodiscard]] inline constexpr auto
is_gc(const char c) {
  return (c >> 1) & 1;
}

inline constexpr auto nibble_size = 4;

[[nodiscard]] inline constexpr auto
encode_nibble(const char c) {
  // ADS: 15 is to keep 4 bits
  // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
  return (c >> 1) & 15;  // N gets separate encoding
}

[[nodiscard]] inline constexpr auto
pct(const double a) {
  return 100.0 * a;
}

inline constexpr auto
add(std::ranges::forward_range auto &a1,
    const std::ranges::forward_range auto &a2) {
  std::ranges::transform(a1, a2, std::begin(a1), std::plus{});
};

inline constexpr auto
add(std::ranges::forward_range auto &a1,
    const std::ranges::forward_range auto &a2, const auto adder) {
  std::ranges::transform(a1, a2, std::begin(a1), adder);
};

template <typename T>
concept has_plus_equals = std::regular<T> && requires(T x, T y) {
  { x += y } -> std::same_as<T &>;
};

inline constexpr auto
add(has_plus_equals auto &a1, const has_plus_equals auto &a2) {
  a1 += a2;
};

inline constexpr auto
vec_add(auto &v1, const auto &v2) {
  v1.resize(std::max(std::size(v1), std::size(v2)));
  add(v1, v2);
};

inline constexpr auto
two_dim_add(auto &v1, const auto &v2) {
  v1.resize(std::max(std::size(v1), std::size(v2)));
  for (auto [a1, a2] : std::views::zip(v1, v2))
    add(a1, a2);
};

inline constexpr auto
two_dim_vec_add(auto &v1, const auto &v2) {
  v1.resize(std::max(std::size(v1), std::size(v2)));
  for (auto [a1, a2] : std::views::zip(v1, v2))
    vec_add(a1, a2);
};

[[nodiscard]] inline constexpr auto
as_frac(const auto a, const auto b) {
  return static_cast<double>(a) / static_cast<double>(b);
}

[[nodiscard]] consteval auto
ipow(const auto b, const auto e) -> std::remove_cvref_t<decltype(b)> {
  return e == 0 ? 1 : (e & 1 ? b : 1) * ipow(b * b, e >> 1);
}

inline constexpr auto
count_nucs(auto seq_itr, const auto seq_end,
           auto &tab) {  // cppcheck-suppress constParameterReference
  auto out_itr = std::begin(tab);
  while (seq_itr != seq_end)
    ++(*out_itr++)[encode(*seq_itr++)];
}

inline constexpr auto
count_ns(auto seq_itr, const auto seq_end,
         auto &tab) {  // cppcheck-suppress constParameterReference
  auto out_itr = std::begin(tab);
  while (seq_itr != seq_end)
    *out_itr++ += (*seq_itr++ == 'N');
}

[[nodiscard]] static inline auto
count_gc(auto seq_itr, const auto seq_end) {
  auto gc = 0;
  while (seq_itr != seq_end) {
    gc += (*seq_itr != 'N' && is_gc(*seq_itr));
    ++seq_itr;
  }
  return gc;
}

[[nodiscard]] static inline auto
count_quals(auto qual_itr, const auto qual_end,
            auto &tab) {  // cppcheck-suppress constParameterReference
  auto out_itr = std::begin(tab);
  auto qual_tot = 0;
  while (qual_itr != qual_end) {
    const auto q = *qual_itr++;
    ++(*out_itr++)[q];
    qual_tot += q;
  }
  return qual_tot;
}

[[nodiscard]] static inline auto
count_quals_rev(auto qual_itr, const auto qual_end,
                auto &tab) {  // cppcheck-suppress constParameterReference
  auto out_itr = std::begin(tab) + std::distance(qual_itr, qual_end);
  auto qual_tot = 0;
  while (qual_itr != qual_end) {
    const auto q = *qual_itr++;
    ++(*(--out_itr))[q];
    qual_tot += q;
  }
  return qual_tot;
}

static inline auto
count_quals_itr(auto qual_itr, const auto qual_end, auto tab_itr) {
  while (qual_itr != qual_end) {
    tab_itr->first += *qual_itr++;
    ++tab_itr->second;
    ++tab_itr;
  }
}

static inline auto
count_quals_itr_rev(auto qual_itr, const auto qual_end, auto tab_itr) {
  tab_itr += std::distance(qual_itr, qual_end);
  while (qual_itr != qual_end) {
    --tab_itr;
    tab_itr->first += *qual_itr++;
    ++tab_itr->second;
  }
}

[[nodiscard]] inline auto
tabular_dot(const auto &a) {
  auto total = static_cast<std::remove_cvref_t<decltype(a)>::value_type>(0);
  for (const auto [i, x] : std::views::enumerate(a))
    total += i * x;
  return total;
}

[[nodiscard]] inline auto
mean_tabular(const auto &a) {
  const auto num = tabular_dot(a);
  const auto denom = std::reduce(std::cbegin(a), std::cend(a));
  return static_cast<double>(num) / static_cast<double>(denom);
}

[[nodiscard]] inline auto
median_tabular(const auto &a) {
  using value_type = std::decay_t<decltype(a)>::value_type;
  std::vector<value_type> cumul(std::size(a), 0);
  std::inclusive_scan(std::cbegin(a), std::cend(a), std::begin(cumul));
  const auto ub = std::ranges::upper_bound(cumul, cumul.back() / 2);
  return static_cast<std::uint32_t>(std::distance(std::begin(cumul), ub));
}

// clang-format off
// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
[[nodiscard]] inline constexpr auto median_val(const auto &q) { return q[0]; }
[[nodiscard]] inline constexpr auto lquart_val(const auto &q) { return q[1]; }
[[nodiscard]] inline constexpr auto uquart_val(const auto &q) { return q[2]; }
[[nodiscard]] inline constexpr auto ldec_val(const auto &q) { return q[3]; }
[[nodiscard]] inline constexpr auto udec_val(const auto &q) { return q[4]; }
// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
// clang-format on

// NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
[[nodiscard]] inline auto
five_quants(const auto &a) -> std::array<std::uint32_t, 5> {
  const auto dist_to_insertion_point = [](const auto &p, const auto x) {
    // get quantile as distance to insertion point in cumulative counts; using
    // upper bound for the case of very few counts
    const auto ub = std::ranges::upper_bound(p, x);
    return static_cast<std::uint32_t>(std::distance(std::begin(p), ub));
  };
  std::vector<std::uint64_t> cumul(std::size(a), 0);
  std::inclusive_scan(std::cbegin(a), std::cend(a), std::begin(cumul));
  return {
    dist_to_insertion_point(cumul, cumul.back() / 2),      // median
    dist_to_insertion_point(cumul, cumul.back() / 4),      // lower quartile
    dist_to_insertion_point(cumul, 3 * cumul.back() / 4),  // upper quartile
    dist_to_insertion_point(cumul, cumul.back() / 10),     // 10th percentile
    dist_to_insertion_point(cumul, 9 * cumul.back() / 10)  // 90th percentile
  };
}
// NOLINTEND(cppcoreguidelines-avoid-magic-numbers)

[[nodiscard]] auto
size_to_units(const std::int64_t s,
              const std::string &suffix = "iB") -> std::string;

[[nodiscard]] static inline auto
get_max_size(const auto &x) {
  assert(!x.empty());
  const auto sz = [](const auto &y) { return std::size(y); };
  return std::ranges::max(std::views::transform(x | std::views::values, sz));
}

#endif  // SRC_FALCO_UTILS_HPP_
