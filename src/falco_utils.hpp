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

#include "falco_file_format.hpp"
#include "quality_score.hpp"

#include <config.h>

#include <htslib/hts.h>
#include <htslib/thread_pool.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <format>
#include <iterator>
#include <numeric>
#include <ranges>
#include <vector>

struct run_mode {
  bool do_tiles{};
  bool do_kmers{};

  // clang-format off
  auto tiles(const bool x) { do_tiles = x; }
  [[nodiscard]] auto tiles() const { return do_tiles; }
  auto kmers(const bool x) { do_kmers = x; }
  [[nodiscard]] auto kmers() const { return do_kmers; }
  [[nodiscard]] auto string() const -> std::string { return {}; }
  // clang-format on
};

// clang-format off
[[nodiscard]] constexpr inline auto tiles(const run_mode &rm) { return rm.tiles(); }
[[nodiscard]] constexpr inline auto kmers(const run_mode &rm) { return rm.kmers(); }
// clang-format on

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
  return (c >> 1) & 15;  // N gets separate encoding
}

inline constexpr auto pct = [](const double a) { return 100.0 * a; };

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

[[nodiscard]] consteval auto
ipow(const auto b, const auto e) -> std::remove_cvref_t<decltype(b)> {
  return e == 0 ? 1 : (e & 1 ? b : 1) * ipow(b * b, e >> 1);
}

static inline auto
count_nucs(auto seq_itr, const auto seq_end,
           auto &tab) {  // cppcheck-suppress constParameterReference
  auto out_itr = std::begin(tab);
  while (seq_itr != seq_end)
    ++(*out_itr++)[encode(*seq_itr++)];
}

static inline auto
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
    ++(*tab_itr++).second;
  }
}

static inline auto
count_quals_itr_rev(auto qual_itr, const auto qual_end, auto tab_itr) {
  tab_itr += std::distance(qual_itr, qual_end);
  while (qual_itr != qual_end) {
    tab_itr->first += *qual_itr++;
    ++(*(--tab_itr)).second;
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

// NOLINTBEGIN (cppcoreguidelines-avoid-magic-numbers)
[[nodiscard]] inline auto
five_quants(const auto &a) -> std::array<std::uint32_t, 5> {
  const auto dlb = [](const auto &p, const auto x) {
    // get quantile as distance to insertion point
    return static_cast<std::uint32_t>(
      std::distance(std::begin(p), std::ranges::lower_bound(p, x)));
  };
  std::vector<std::uint64_t> p(std::size(a), 0);
  std::inclusive_scan(std::cbegin(a), std::cend(a), std::begin(p));
  return {
    dlb(p, p.back() / 2),      // median
    dlb(p, p.back() / 4),      // lower quartile
    dlb(p, 3 * p.back() / 4),  // upper quartile
    dlb(p, p.back() / 10),     // 10th percentile
    dlb(p, 9 * p.back() / 10)  // 90th percentile
  };
}
// NOLINTEND (cppcoreguidelines-avoid-magic-numbers)

[[nodiscard]] inline auto
get_grade(const auto &cutoffs, const auto c) {
  const auto a = std::pair{c, std::string{}};
  const auto b = std::ranges::lower_bound(cutoffs, a);
  if (b == std::cend(cutoffs))
    throw std::runtime_error("error in identifying grade");
  return b->second;
}

struct falco_thread_pool {
  htsThreadPool t{};
  explicit falco_thread_pool(const std::uint32_t n_threads) :
    t{hts_tpool_init(std::max(1u, n_threads)), 0} {
    if (t.pool == nullptr)
      throw std::runtime_error("failed to construct thread pool");
  }
  ~falco_thread_pool() { hts_tpool_destroy(t.pool); }
};

static inline constexpr auto
falco_reverse(auto first, auto last) {
  for (--last; first < last; ++first, --last)
    std::iter_swap(first, last);
}

static inline constexpr auto
falco_complement(auto first, auto last) {
  const auto complement = [](auto &a) {
    a = "TNGNNNCNNNNNNNNNNNNANNNNNN"[a - 'A'];
  };
  while (first != last)
    complement(*first++);
}

static inline constexpr auto
revcomp(auto first, auto last) {
  if (first == last) [[unlikely]]
    return;
  falco_reverse(first, last);
  falco_complement(first, last);
}

#endif  // SRC_FALCO_UTILS_HPP_
