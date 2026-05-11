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

[[nodiscard]] inline auto
format_read_lengths(const auto &lengths, const auto max_read_len) {
  static constexpr auto start_tag = ">>Sequence Length Distribution\t{}\n";
  static constexpr auto header = "#Length Count\n";
  const bool has_empty_reads = std::size(lengths) > 0 && lengths[0] > 0;
  const auto n_lengths =
    std::ranges::count_if(lengths, [](const auto x) { return x > 0; });
  const auto grade = has_empty_reads ? "fail" : n_lengths > 1 ? "warn" : "pass";
  auto r = std::format(start_tag, grade);
  r += header;
  for (auto i = 0u; i <= max_read_len; ++i)
    if (lengths[i] > 0)
      r += std::format("{}\t{}\n", i, lengths[i]);
  r += end_module_tag;
  return r;
}

// calculate deviation of a hist from a "normal" with same mode and sd
[[nodiscard]] static inline auto
sum_deviation_from_normal(const auto &gc) {
  const auto n_bins = std::size(gc);
  const auto gc_beg = std::cbegin(gc);
  const auto gc_end = std::cend(gc);
  const auto total_count = std::reduce(gc_beg, gc_end);
  assert(n_bins > 1 && total_count > 1);

  // get mode
  const auto mode_itr = std::ranges::max_element(gc);
  const std::uint64_t mode_pos = std::distance(gc_beg, mode_itr);
  const auto mode_val = *mode_itr;

  // ADS: in case mode is not sharp average nearby values (not clear on why)
  const auto gt_cut = [&](const double x) { return x < (0.90 * mode_val); };
  const auto right_itr = std::find_if(mode_itr, std::cend(gc), gt_cut);
  const auto left_itr =
    std::find_if(std::reverse_iterator(mode_itr), std::crend(gc), gt_cut);
  const auto n = std::distance(std::reverse_iterator(right_itr), left_itr);
  const auto mode = mode_pos + (n - 1) / 2.0;

  // theoretical distribution
  const auto cntr_sq = [m = mode](const auto z) { return (z - m) * (z - m); };
  const auto sd_term = [&](const auto &x) {
    return cntr_sq(std::get<0>(x)) * std::get<1>(x);
  };
  const auto id_gc = std::views::enumerate(gc) | std::views::transform(sd_term);
  const auto sd = std::sqrt(std::reduce(std::cbegin(id_gc), std::cend(id_gc))) /
                  (total_count - 1);

  const auto norm = [&](const auto i) {
    return std::exp(-cntr_sq(i) / (2.0 * sd * sd));
  };
  auto th = std::views::iota(0u, n_bins) | std::views::transform(norm) |
            std::ranges::to<std::vector>();
  const auto tot = std::reduce(std::cbegin(th), std::cend(th));
  std::ranges::transform(th, std::begin(th),
                         [&](const auto x) { return x * total_count / tot; });
  const auto diff = [](const auto a, const auto b) { return std::fabs(a - b); };
  const auto r = std::transform_reduce(gc_beg, gc_end, std::cbegin(th), 0.0,
                                       std::plus{}, diff);

  return pct(as_frac(r, total_count));
}

[[nodiscard]] inline auto
format_gc_content(const auto &gc_content, const auto max_read_len) {
  static constexpr auto start_tag = ">>Per sequence GC content\t{}\n";
  static constexpr auto header = "#GC\tContent\tCount\n";
  static constexpr auto grade_cutoffs = std::array{
    std::pair{15.0, "pass"},
    std::pair{30.0, "warn"},
    std::pair{std::numeric_limits<double>::max(), "error"},
  };
  const auto deviation = sum_deviation_from_normal(gc_content);
  auto r = std::format(start_tag, get_grade(grade_cutoffs, deviation));
  r += header;
  for (auto i = 0u; i < max_read_len; ++i)
    // cppcheck-suppress useStlAlgorithm
    r += std::format("{}\t{}\n", i + 1, gc_content[i]);
  r += end_module_tag;

  return r;
}

[[nodiscard]] inline auto
format_base_composition(const auto &nucs, const auto max_read_len) {
  static constexpr auto start_tag = ">>Per base sequence content\t{}\n";
  static constexpr auto header = "#Base\tG\tA\tT\tC\n";
  static constexpr auto base_permutation = {3, 0, 2, 1};
  auto r = std::format(start_tag, "fail");
  r += header;
  for (auto i = 0u; i < max_read_len; ++i) {
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

[[nodiscard]] inline auto
format_n_counts(const auto &n_counts, const auto &nucs,
                const auto max_read_len) {
  static constexpr auto start_tag = "Per base N content\t{}\n";
  static constexpr auto header = "#Base\tN-Count\n";
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.05, "pass"},
    std::pair{0.20, "warn"},
    std::pair{1.00, "error"},
  };
  const auto max_idx =
    std::distance(std::cbegin(n_counts), std::ranges::max_element(n_counts));
  const auto max_n =
    as_frac(n_counts[max_idx],
            std::reduce(std::cbegin(nucs[max_idx]), std::cend(nucs[max_idx])));
  auto r = std::format(start_tag, get_grade(grade_cutoffs, max_n));
  r += header;
  for (auto i = 0u; i < max_read_len; ++i) {
    const auto tot = std::reduce(std::cbegin(nucs[i]), std::cend(nucs[i]));
    r += std::format("{}\t{:.6g}\n", i + 1, pct(as_frac(n_counts[i], tot)));
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] inline auto
format_qual_by_read(const auto &qual_by_read, const auto qual_offset) {
  static constexpr auto start_tag = ">>Per sequence quality scores\t{}\n";
  static constexpr auto header = "#Quality\tCount\n";
  static constexpr auto grade_cutoffs = std::array{
    std::pair{20, "error"},
    std::pair{27, "warn"},
    std::pair{falco::max_qual_val, "pass"},
  };
  // get mode
  const auto max_qual_val =
    std::distance(std::cbegin(qual_by_read),
                  std::ranges::max_element(qual_by_read)) -
    qual_offset;
  auto r = std::format(start_tag, get_grade(grade_cutoffs, max_qual_val));
  r += header;
  // output starting at qual_offset; that's where they are relevant
  for (const auto q : std::views::iota(qual_offset, falco::max_qual_val))
    if (qual_by_read[q] > 0)
      r += std::format("{}\t{}\n", q - qual_offset, qual_by_read[q]);
  r += end_module_tag;
  return r;
}

[[nodiscard]] inline auto
get_grade_qual_by_pos(const auto &qual, const auto max_read_len) {
  static constexpr auto get_lower_quartile = [](const auto &quantiles) {
    return quantiles[1];  // NOLINT (*-avoid-magic-numbers)
  };
  static constexpr auto get_median = [](const auto &quantiles) {
    return quantiles[2];  // NOLINT (*-avoid-magic-numbers)
  };
  // ADS: not sure how this is used
  [[maybe_unused]] static constexpr auto grade_cutoffs_lower = std::array{
    std::pair{0.05, "error"},
    std::pair{0.10, "warn"},
    std::pair{1.00, "pass"},
  };
  static constexpr auto grade_cutoffs_median = std::array{
    std::pair{0.20, "error"},
    std::pair{0.25, "warn"},
    std::pair{1.00, "pass"},
  };
  auto min_qual_median = 1.0;
  for (const auto &q : qual | std::views::take(max_read_len)) {
    const auto quantiles = five_quants(q);
    // ADS: lower quartile to ensure enough data?
    const auto median =
      get_lower_quartile(quantiles) > 0.0 ? get_median(quantiles) : 1.0;
    min_qual_median = std::min(min_qual_median, median);
  }
  return get_grade(grade_cutoffs_median, min_qual_median);
}

[[nodiscard]] inline auto
format_qual_by_pos(const auto &qual, const auto max_read_len,
                   const auto qual_offset) {
  static constexpr auto start_tag = ">>Per base sequence quality\t{}\n";
  static constexpr auto header = "#Base\tMean\tMedian\tLower Quartile\tUpper "
                                 "Quartile\t10th Percentile\t90th Percentile\n";
  const auto tab_sep = [](const auto &a) {
    auto r = std::string{};
    for (const auto &value : a | std::views::take(std::size(a) - 1))
      // cppcheck-suppress useStlAlgorithm
      r += std::format("{}\t", value);
    return r + std::format("{}", a.back());
  };
  const auto sub_from_each = [](auto a, const auto b) {
    std::ranges::for_each(a, [&](auto &x) { x -= b; });
    return a;
  };
  auto r = std::format(start_tag, get_grade_qual_by_pos(qual, max_read_len));
  r += header;
  for (auto i = 0u; i < max_read_len; ++i) {
    const auto q = sub_from_each(five_quants(qual[i]), qual_offset);
    r += std::format("{}\t{:.6g}\t{}\n", i + 1,
                     mean_tabular(qual[i]) - qual_offset, tab_sep(q));
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] inline auto
format_basic_stats(const auto &filename, const auto n_reads,
                   const auto max_read_len, const auto total_gc,
                   const auto total_nucs, const auto &encoding) {
  static constexpr auto grade = "pass";  // always a pass
  static constexpr auto start_tag = "##Falco {}\n"
                                    ">>Basic Statistics\t{}\n";
  static constexpr auto header = "#Measure\tValue\n";
  [[maybe_unused]] const auto with_suffix = [&](const auto x) {
    if (x > 1'000'000'000)
      return std::format("{:.1f} {}bp", as_frac(x, 1'000'000'000), "G");
    if (x > 1'000'000)
      return std::format("{:.1f} {}bp", as_frac(x, 1'000'000), "M");
    return std::format("{:.1f} {}bp", as_frac(x, 1'000), "K");
  };
  auto r = std::format(start_tag, VERSION, grade);
  r += header;
  r += std::format("Filename\t{}\n", filename);
  const auto [_, file_type] = get_file_format(filename);
  r += std::format("File type\t{}\n", file_type);
  r += std::format("Encoding\t{}\n", encoding);
  r += std::format("Total Sequences\t{}\n", n_reads);
  r += std::format("Total Bases\t{}\n", total_nucs);
  r += std::format("Sequences flagged as poor quality {}\n", 0);
  r += std::format("Sequence length\t{}\n", max_read_len);
  r += std::format("%GC\t{:.1f}\n", pct(as_frac(total_gc, total_nucs)));
  r += end_module_tag;
  return r;
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

#endif  // SRC_FALCO_UTILS_HPP_
