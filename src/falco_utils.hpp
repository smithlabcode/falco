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

#include "quality_score.hpp"

#include <htslib/hts.h>
#include <htslib/thread_pool.h>

#include <config.h>

#include <algorithm>
#include <array>
#include <format>
#include <iterator>
#include <numeric>
#include <ranges>
#include <vector>

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

// NOLINTBEGIN (cppcoreguidelines-avoid-magic-numbers)
[[nodiscard]] static inline auto
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

[[nodiscard]] static inline auto
format_read_lengths(const auto &lengths, const auto max_read_len) {
  static constexpr auto start_tag = ">>Sequence Length Distribution\t{}\n";
  static constexpr auto header = "#Length Count\n";
  auto r = std::format(start_tag, "pass");
  r += header;
  for (auto i = 0; i <= max_read_len; ++i)
    if (lengths[i] > 0)
      r += std::format("{}\t{}\n", i, lengths[i]);
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
format_gc_content(const auto &gc_content, const auto max_read_len) {
  static constexpr auto start_tag = ">>Per sequence GC content\t{}\n";
  static constexpr auto header = "#GC\tContent\tCount\n";
  auto r = std::format(start_tag, "pass");
  r += header;
  for (auto i = 0; i < max_read_len; ++i)
    r += std::format("{}\t{}\n", i + 1, gc_content[i]);
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
format_base_composition(const auto &nucs, const auto max_read_len) {
  static constexpr auto start_tag = ">>Per base sequence content\t{}\n";
  static constexpr auto header = "#Base\tG\tA\tT\tC\n";
  static constexpr auto base_permutation = {3, 0, 2, 1};
  auto r = std::format(start_tag, "fail");
  r += header;
  for (auto i = 0; i < max_read_len; ++i) {
    r += std::format("{}", i + 1);
    const auto tot = std::reduce(std::cbegin(nucs[i]), std::cend(nucs[i]));
    for (const auto j : base_permutation)
      // cppcheck-suppress useStlAlgorithm
      r += std::format("\t{:2.4f}", pct(as_frac(nucs[i][j], tot)));
    r += '\n';
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
format_n_counts(const auto &n_counts, const auto &nucs,
                const auto max_read_len) {
  static constexpr auto start_tag = "Per base N content\t{}\n";
  static constexpr auto header = "#Base\tN-Count\n";
  auto r = std::format(start_tag, "pass");
  r += header;
  for (auto i = 0; i < max_read_len; ++i) {
    const auto tot = std::reduce(std::cbegin(nucs[i]), std::cend(nucs[i]));
    r += std::format("{}\t{:.6g}\n", i + 1, pct(as_frac(n_counts[i], tot)));
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
format_qual_by_read(const auto &qual, const auto max_read_len) {
  static constexpr auto start_tag = ">>Per sequence quality scores\t{}\n";
  static constexpr auto header = "#Quality\tCount\n";
  const auto qual_itr = std::cbegin(qual);
  const auto qual_tot =
    std::reduce(qual_itr + falco::min_qual_val, qual_itr + falco::max_qual_val);
  auto r = std::format(start_tag, "pass");
  r += header;
  for (auto i = falco::min_qual_val; i < falco::max_qual_val; ++i)
    r += std::format("{}\t{:.6g}\n", i, as_frac(qual[i], qual_tot));
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
format_qual_by_pos(const auto &qual, const auto max_read_len) {
  static constexpr auto start_tag = ">>Per base sequence quality\t{}\n";
  static constexpr auto header = "#Base\tMean\tMedian\tLower Quartile\tUpper "
                                 "Quartile\t10th Percentile\t90th Percentile\n";
  const auto tab_sep = [](const auto &a) {
    auto r = std::string{};
    for (const auto &value : a | std::views::take(std::size(a) - 1))
      r += std::format("{}\t", value);
    return r + std::format("{}", a.back());
  };
  const auto sub_from_each = [](auto a, const auto b) {
    std::ranges::for_each(a, [&](auto &x) { x -= b; });
    return a;
  };
  auto r = std::format(start_tag, "pass");
  r += header;
  for (auto i = 0; i < max_read_len; ++i) {
    const auto q = sub_from_each(five_quants(qual[i]), falco::min_qual_val);
    r += std::format("{}\t{:.6g}\t{}\n", i + 1,
                     mean_tabular(qual[i]) - falco::min_qual_val, tab_sep(q));
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
format_basic_stats(/*const auto &filename,*/ const auto n_reads,
                   const auto max_read_len, const auto total_gc,
                   const auto total_nucs) {
  static constexpr auto start_tag = "##Falco {}\n>>Basic Statistics\n";
  static constexpr auto header = "#Measure\tValue\n";
  auto r = std::format("##Falco {}\n", VERSION);
  r += header;
  r += std::format("Filename\t{}\n", "asdf");
  r += std::format("File type\t{}\n", "Conventional base calls");
  r += std::format("Encoding\t{}\n", "Sanger / Illumina 1.9");
  r += std::format("Total Sequences\t{}\n", n_reads);
  r += std::format("Sequences flagged as poor quality {}\n", 0);
  r += std::format("Sequence length\t{}\n", max_read_len);
  r += std::format("%GC\t{:.1f}\n", pct(as_frac(total_gc, total_nucs)));
  r += end_module_tag;
  return r;
}

struct hts_thread_pool_wrapper {
  htsThreadPool t{};
  explicit hts_thread_pool_wrapper(const std::uint32_t n_threads) :
    t{hts_tpool_init(std::max(1u, n_threads)), 0} {
    if (t.pool == nullptr)
      throw std::runtime_error("failed to construct thread pool");
  }
  ~hts_thread_pool_wrapper() { hts_tpool_destroy(t.pool); }
};

#endif  // SRC_FALCO_UTILS_HPP_
