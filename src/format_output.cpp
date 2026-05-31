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

#include "format_output.hpp"
#include "falco_utils.hpp"
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
#include <limits>  // IWYU pragma: keep
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

[[nodiscard]] auto
get_grade_read_lengths(const std::vector<std::uint64_t> &lengths)
  -> std::string {
  const bool has_empty_reads = std::size(lengths) > 0 && lengths[0] > 0;
  const auto n_lengths =
    std::ranges::count_if(lengths, [](const auto x) { return x > 0; });
  const auto grade = has_empty_reads ? "fail" : n_lengths > 1 ? "warn" : "pass";
  return grade;
}

[[nodiscard]] auto
format_read_lengths(const std::vector<std::uint64_t> &lengths) -> std::string {
  static constexpr auto start_tag = ">>Sequence Length Distribution\t{}\n";
  static constexpr auto header = "#Length\t"
                                 "Count\n";
  const auto eq0 = [](const auto &x) { return std::get<1>(x) == 0; };
  auto r = std::format(start_tag, get_grade_read_lengths(lengths));
  r += header;
  auto to_report = std::views::enumerate(lengths) | std::views::drop_while(eq0);
  for (const auto [idx, len] : to_report)
    r += std::format("{}\t{}\n", idx, len);
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
sum_deviation_from_normal(const auto &gc) {
  static constexpr auto mode_width = 0.90;  // ADS: where does this come from?
  // calculate deviation of a hist from a "normal" with same mode and sd
  const auto n_bins = std::size(gc);
  const auto gc_beg = std::cbegin(gc);
  const auto gc_end = std::cend(gc);
  const auto total_count = std::reduce(gc_beg, gc_end);
  if (total_count <= 1)  // we will be dividing by (total_count - 1)
    return 0.0;

  // get mode
  const auto mode_itr = std::ranges::max_element(gc);
  const std::uint64_t mode_pos = std::distance(gc_beg, mode_itr);
  const auto mode_val = *mode_itr;

  // ADS: in case mode is not sharp average nearby values (not clear on why)
  const auto gt_cut = [&](const double x) { return x < mode_width * mode_val; };
  const auto right_itr = std::find_if(mode_itr, std::cend(gc), gt_cut);
  const auto left_itr =
    std::find_if(std::reverse_iterator(mode_itr), std::crend(gc), gt_cut);
  const auto n = std::distance(std::reverse_iterator(right_itr), left_itr);
  const double mode = mode_pos + (n - 1) / 2.0;

  // theoretical distribution
  const auto cntr_sq = [m = mode](const auto v) { return (v - m) * (v - m); };
  const auto sd_term = [&](const auto &x) {
    return cntr_sq(std::get<0>(x)) * std::get<1>(x);
  };
  const auto id_gc = std::views::enumerate(gc) | std::views::transform(sd_term);
  const auto sd = std::sqrt(std::reduce(std::cbegin(id_gc), std::cend(id_gc)) /
                            (total_count - 1));
  const auto to_normal = [&](const auto val) {
    // NOLINTNEXTLINE (cppcoreguidelines-avoid-magic-numbers)
    return std::exp(-cntr_sq(val) / (2.0 * sd * sd));
  };
  auto theor = std::views::iota(0u, n_bins) | std::views::transform(to_normal) |
               std::ranges::to<std::vector>();
  const auto denom = std::reduce(std::cbegin(theor), std::cend(theor));
  std::ranges::transform(theor, std::begin(theor),
                         [&](const auto x) { return x * total_count / denom; });
  const auto diff = [](const auto a, const auto b) { return std::fabs(a - b); };
  const auto r = std::transform_reduce(gc_beg, gc_end, std::cbegin(theor), 0.0,
                                       std::plus{}, diff);
  return pct(as_frac(r, total_count));
}

[[nodiscard]] auto
get_grade_gc_content(const falco::gc_content_array &gc_content) -> std::string {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{15.0, "pass"},
    std::pair{30.0, "warn"},
    std::pair{std::numeric_limits<double>::max(), "fail"},
  };
  const auto deviation = sum_deviation_from_normal(gc_content);
  return get_grade(grade_cutoffs, deviation);
}

[[nodiscard]] auto
format_gc_content(const falco::gc_content_array &gc_content) -> std::string {
  static constexpr auto start_tag = ">>Per sequence GC content\t{}\n";
  static constexpr auto header = "#GC Content\t"
                                 "Count\n";
  auto r = std::format(start_tag, get_grade_gc_content(gc_content));
  r += header;
  for (const auto [idx, gc] : std::views::enumerate(gc_content))
    r += std::format("{}\t{}\n", idx, gc);
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
get_grade_base_composition(const std::vector<falco::nuc_array> &nucs)
  -> std::string {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{10.0, "pass"},
    std::pair{20.0, "warn"},
    std::pair{std::numeric_limits<double>::max(), "fail"},
  };
  const auto compl_diff = [](const auto &by_pos) {
    // (A,C,G,T)=(0,1,3,2)
    const auto tot = std::reduce(std::cbegin(by_pos), std::cend(by_pos));
    const auto delta = [tot](const auto a, const auto b) {
      return pct(as_frac(a, tot)) - pct(as_frac(b, tot));
    };
    return std::max(std::fabs(delta(by_pos[0], by_pos[2])),
                    std::fabs(delta(by_pos[1], by_pos[3])));
  };
  const auto max_diff =
    std::ranges::max(nucs | std::views::transform(compl_diff));
  return get_grade(grade_cutoffs, max_diff);
}

[[nodiscard]] auto
format_base_composition(const std::vector<falco::nuc_array> &nucs)
  -> std::string {
  static constexpr auto start_tag = ">>Per base sequence content\t{}\n";
  static constexpr auto header = "#Base\t"
                                 "G\tA\tT\tC\n";
  static constexpr auto base_permutation = {3, 0, 2, 1};
  auto r = std::format(start_tag, get_grade_base_composition(nucs));
  r += header;
  for (auto i = 0u; i < std::size(nucs); ++i) {
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

[[nodiscard]] auto
get_grade_n_counts(const std::vector<std::uint64_t> &n_counts,
                   const std::vector<falco::nuc_array> &nucs) -> std::string {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.05, "pass"},
    std::pair{0.20, "warn"},
    std::pair{1.00, "fail"},
  };
  const auto max_idx =
    std::distance(std::cbegin(n_counts), std::ranges::max_element(n_counts));
  const auto total =
    std::reduce(std::cbegin(nucs[max_idx]), std::cend(nucs[max_idx]));
  const auto max_n_frac = as_frac(n_counts[max_idx], total);
  return get_grade(grade_cutoffs, max_n_frac);
}

[[nodiscard]] auto
format_n_counts(const std::vector<std::uint64_t> &n_counts,
                const std::vector<falco::nuc_array> &nucs) -> std::string {
  static constexpr auto start_tag = ">>Per base N content\t{}\n";
  static constexpr auto header = "#Base\t"
                                 "N-Count\n";
  auto r = std::format(start_tag, get_grade_n_counts(n_counts, nucs));
  r += header;
  for (auto i = 0u; i < std::size(n_counts); ++i) {
    const auto tot = std::reduce(std::cbegin(nucs[i]), std::cend(nucs[i]));
    r += std::format("{}\t{:.6g}\n", i + 1, pct(as_frac(n_counts[i], tot)));
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
get_grade_qual_by_read(const falco::qual_array &qual_by_read) -> std::string {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{20, "fail"},
    std::pair{27, "warn"},
    std::pair{falco::max_qual_val, "pass"},
  };
  // get mode for grade
  const auto q_beg = std::cbegin(qual_by_read);
  const auto max_itr = std::ranges::max_element(qual_by_read);
  const auto qual_val_mode = std::distance(q_beg, max_itr);
  return get_grade(grade_cutoffs, qual_val_mode);
}

[[nodiscard]] auto
format_qual_by_read(const falco::qual_array &qual_by_read) -> std::string {
  static constexpr auto start_tag = ">>Per sequence quality scores\t{}\n";
  static constexpr auto header = "#Quality\t"
                                 "Count\n";
  auto r = std::format(start_tag, get_grade_qual_by_read(qual_by_read));
  r += header;

  // output quality values between first non-zero and last zero
  const auto gt0 = [&](const auto x) { return x > 0; };
  const auto first_obs_itr = std::ranges::find_if(qual_by_read, gt0);
  const auto q_beg = std::cbegin(qual_by_read);
  const std::int64_t first_obs = std::distance(q_beg, first_obs_itr);
  const auto last_obs_subrange = std::ranges::find_last_if(qual_by_read, gt0);
  const auto trailing_zeros = std::size(last_obs_subrange);
  const std::int64_t last_obs = std::size(qual_by_read) - trailing_zeros;
  assert(first_obs >= 0 && last_obs <= falco::max_qual_val);
  for (const auto q : std::views::iota(first_obs, last_obs + 1))
    // cppcheck-suppress useStlAlgorithm
    r += std::format("{}\t{}\n", q, qual_by_read[q]);
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
get_grade_qual_by_pos(const std::vector<falco::qual_array> &qual)
  -> std::string {
  // NOLINTBEGIN (cppcoreguidelines-avoid-magic-numbers)
  static constexpr auto lower_quartile = [](const auto &q) { return q[1]; };
  static constexpr auto get_median = [](const auto &q) { return q[2]; };
  // NOLINTEND (cppcoreguidelines-avoid-magic-numbers)
  // ADS: not sure how this is used
  [[maybe_unused]] static constexpr auto grade_cutoffs_lower = std::array{
    std::pair{0.05, "fail"},
    std::pair{0.10, "warn"},
    std::pair{1.00, "pass"},
  };
  static constexpr auto grade_cutoffs_median = std::array{
    std::pair{0.20, "fail"},
    std::pair{0.25, "warn"},
    std::pair{1.00, "pass"},
  };
  auto min_qual_median = 1.0;
  for (const auto &q : qual) {
    const auto quantiles = five_quants(q);
    // ADS: is lower quartile to ensure enough data?
    const auto median =
      lower_quartile(quantiles) > 0.0 ? get_median(quantiles) : 1.0;
    min_qual_median = std::min(min_qual_median, median);
  }
  return get_grade(grade_cutoffs_median, min_qual_median);
}

[[nodiscard]] auto
format_qual_by_pos(const std::vector<falco::qual_array> &qual) -> std::string {
  static constexpr auto digits{std::numeric_limits<double>::digits10};
  static constexpr auto start_tag = ">>Per base sequence quality\t{}\n";
  static constexpr auto header = "#Base\t"
                                 "Mean\t"
                                 "Median\t"
                                 "Lower Quartile\t"
                                 "Upper Quartile\t"
                                 "10th Percentile\t"
                                 "90th Percentile\n";
  const auto tab_sep = [](const auto a) {
    const auto with_tab = [](const auto x) { return std::format("\t{}", x); };
    return std::views::transform(a, with_tab) | std::views::join |
           std::ranges::to<std::string>();
  };
  auto r = std::format(start_tag, get_grade_qual_by_pos(qual));
  r += header;
  for (const auto [idx, q] : std::views::enumerate(qual))
    r += std::format("{}\t{:.{}g}{}\n", idx + 1, mean_tabular(q), digits,
                     tab_sep(five_quants(q)));
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
format_basic_stats(const std::string &filename, const std::uint64_t n_reads,
                   const std::uint64_t min_read_len,
                   const std::uint64_t max_read_len,
                   const std::uint64_t total_gc, const std::uint64_t total_nucs,
                   const std::string &encoding) -> std::string {
  static constexpr auto grade = "pass";  // always a pass
  static constexpr auto start_tag = "##Falco {}\n"
                                    ">>Basic Statistics\t{}\n";
  static constexpr auto header = "#Measure\t"
                                 "Value\n";
  [[maybe_unused]] const auto with_suffix = [&](const auto x) {
    // NOLINTBEGIN (cppcoreguidelines-avoid-magic-numbers)
    if (x > 1'000'000'000)
      return std::format("{:.1f} {}bp", as_frac(x, 1'000'000'000), "G");
    if (x > 1'000'000)
      return std::format("{:.1f} {}bp", as_frac(x, 1'000'000), "M");
    return std::format("{:.1f} {}bp", as_frac(x, 1'000), "K");
    // NOLINTEND (cppcoreguidelines-avoid-magic-numbers)
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
  r += std::format("Sequence length\t{}\n",
                   min_read_len == max_read_len
                     ? std::format("{}", max_read_len)
                     : std::format("{}-{}", min_read_len, max_read_len));
  r += std::format("%GC\t{:.1f}\n", pct(as_frac(total_gc, total_nucs)));
  r += end_module_tag;
  return r;
}
