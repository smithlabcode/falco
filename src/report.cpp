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

#include "report.hpp"

#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "kmer_counter.hpp"  // IWYU pragma: keep
#include "quality_score.hpp"
#include "tile_processor.hpp"

#define FMT_HEADER_ONLY
#include "fmt/base.h"
#include "fmt/format.h"

#include <config.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <format>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

[[nodiscard]] auto
sequence_length_report(const std::vector<std::uint64_t> &lengths,
                       const file_grades &grades) -> std::string {
  static constexpr auto label = "sequence_length";
  static constexpr auto start_tag = ">>Sequence Length Distribution\t{}\n";
  static constexpr auto header = "#Length\t"
                                 "Count\n";
  const auto eq0 = [](const auto &x) { return std::get<1>(x) == 0; };
  auto r = std::format(start_tag, grades.grade(label));
  r += header;
  auto to_report = std::views::enumerate(lengths) | std::views::drop_while(eq0);
  for (const auto [idx, len] : to_report)
    r += std::format("{}\t{}\n", idx, len);
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
gc_sequence_report(const falco::gc_content_array &gc_content,
                   const file_grades &grades) -> std::string {
  static constexpr auto label = "gc_sequence";
  static constexpr auto start_tag = ">>Per sequence GC content\t{}\n";
  static constexpr auto header = "#GC Content\t"
                                 "Count\n";
  auto r = std::format(start_tag, grades.grade(label));
  r += header;
  for (const auto [idx, gc] : std::views::enumerate(gc_content))
    r += std::format("{}\t{}\n", idx, gc);
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
sequence_report(const std::vector<falco::nuc_array> &nucs,
                const std::vector<base_group_t> &groups,
                const file_grades &grades) -> std::string {
  static constexpr auto label = "sequence";
  static constexpr auto start_tag = ">>Per base sequence content\t{}\n";
  static constexpr auto header = "#Base\t"
                                 "G\tA\tT\tC\n";
  static constexpr auto base_permutation = {3, 0, 2, 1};
  auto r = std::format(start_tag, grades.grade(label));
  r += header;
  for (auto i = 0u; i < std::size(nucs); ++i) {
    r += make_group_tag(groups[i]);
    const auto tot = std::reduce(std::cbegin(nucs[i]), std::cend(nucs[i]));
    for (const auto j : base_permutation)
      // cppcheck-suppress useStlAlgorithm
      // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
      r += std::format("\t{:2.4f}", pct(as_frac(nucs[i][j], tot)));
    r += '\n';
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
n_content_report(const std::vector<std::uint64_t> &n_counts,
                 const std::vector<falco::nuc_array> &nucs,
                 const std::vector<base_group_t> &groups,
                 const file_grades &grades) -> std::string {
  static constexpr auto label = "n_content";
  static constexpr auto start_tag = ">>Per base N content\t{}\n";
  static constexpr auto header = "#Base\t"
                                 "N-Count\n";
  auto r = std::format(start_tag, grades.grade(label));
  r += header;
  for (auto i = 0u; i < std::size(n_counts); ++i) {
    const auto total_non_n =
      std::reduce(std::cbegin(nucs[i]), std::cend(nucs[i]));
    r += std::format("{}\t{:.6g}\n", make_group_tag(groups[i]),
                     pct(as_frac(n_counts[i], n_counts[i] + total_non_n)));
  }
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
quality_sequence_report(const falco::qual_array &qual_by_read,
                        const file_grades &grades) -> std::string {
  static constexpr auto label = "quality_sequence";
  static constexpr auto start_tag = ">>Per sequence quality scores\t{}\n";
  static constexpr auto header = "#Quality\t"
                                 "Count\n";
  auto r = std::format(start_tag, grades.grade(label));
  r += header;

  // output quality values between first non-zero and last zero
  const auto gt0 = [&](const auto x) { return x > 0; };
  const auto first_obs_itr = std::ranges::find_if(qual_by_read, gt0);
  const auto q_beg = std::cbegin(qual_by_read);
  const std::int64_t first_obs = std::distance(q_beg, first_obs_itr);
  const auto last_obs_subrange = std::ranges::find_last_if(qual_by_read, gt0);
  const auto trailing_zeros = std::ssize(last_obs_subrange);
  const std::int64_t last_obs = std::ssize(qual_by_read) - trailing_zeros;
  assert(first_obs >= 0 && last_obs <= falco::max_qual_val);
  for (const auto q : std::views::iota(first_obs, last_obs + 1))
    // cppcheck-suppress useStlAlgorithm
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
    r += std::format("{}\t{}\n", q, qual_by_read[q]);
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
quality_base_report(const std::vector<falco::qual_array> &qual,
                    const std::vector<base_group_t> &groups,
                    const file_grades &grades) -> std::string {
  static constexpr auto digits{std::numeric_limits<double>::digits10};
  static constexpr auto label = "quality_base";
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
  auto r = std::format(start_tag, grades.grade(label));
  r += header;
  for (const auto [idx, q] : std::views::enumerate(qual))
    r += std::format("{}\t{:.{}g}{}\n", make_group_tag(groups[idx]),
                     mean_tabular(q), digits, tab_sep(five_quants(q)));
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
basic_stats_report(const file_info &info, const std::uint64_t n_reads,
                   const std::uint64_t min_read_len,
                   const std::uint64_t max_read_len,
                   const std::uint64_t median_read_len,
                   const std::uint64_t total_gc, const std::uint64_t total_nucs,
                   const file_grades &grades) -> std::string {
  static constexpr auto label = "basic_stats";
  static constexpr auto start_tag = "##Falco {}\n"
                                    ">>Basic Statistics\t{}\n";
  static constexpr auto header = "#Measure\t"
                                 "Value\n";
  [[maybe_unused]] const auto with_suffix = [&](const auto x) {
    // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
    if (x > 1'000'000'000)
      return std::format("{:.1f} {}bp", as_frac(x, 1'000'000'000), "G");
    if (x > 1'000'000)
      return std::format("{:.1f} {}bp", as_frac(x, 1'000'000), "M");
    return std::format("{:.1f} {}bp", as_frac(x, 1'000), "K");
    // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
  };
  auto r = std::format(start_tag, VERSION, grades.grade(label));
  r += header;
  r += std::format("Filename\t{}\n", info.name);
  r += std::format("File type\t{}\n", info.description);
  r += std::format("Encoding\t{}\n", to_string(info.encoding));
  r += std::format("Total Sequences\t{}\n", n_reads);
  r += std::format("Total Bases\t{}\n", total_nucs);
  r += std::format("Sequence length\t{}\n",
                   min_read_len == max_read_len
                     ? std::format("{}", max_read_len)
                     : std::format("{}-{}", min_read_len, max_read_len));
  r += std::format("Mean Length\t{}\n", as_frac(total_nucs, n_reads));
  r += std::format("Median Length\t{}\n", median_read_len);
  r += std::format("%GC\t{:.1f}\n", pct(as_frac(total_gc, total_nucs)));
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
tile_report(const tile_processor::tiles_centered_t &centered,
            const std::vector<base_group_t> &groups,
            const file_grades &grades) -> std::string {
  static constexpr auto label = "tile";
  static constexpr auto max_precision{std::numeric_limits<double>::digits10};
  static constexpr auto start_tag = ">>Per tile sequence quality\t{}\n";
  static constexpr auto header = "#Tile\t"
                                 "Base\t"
                                 "Mean\n";
  if (centered.empty())
    return {};  // ADS: if this fun is called for files w/o tile info
  auto r = std::format(start_tag, grades.grade(label));
  r += header;
  for (const auto &[tile_id, q] : centered)
    for (auto j = 0u; j < std::size(q); ++j)
      r += std::format("{}\t{}\t{:.{}f}\n", tile_id, make_group_tag(groups[j]),
                       q[j], max_precision);
  return r + end_module_tag;
}

[[nodiscard]] auto
kmer_report(const std::vector<kmer_result> &results,
            const file_grades &grades) -> std::string {
  static constexpr auto label = "kmer";
  static constexpr auto start_tag = ">>Kmer Content\t{}\n";
  static constexpr auto header = "#Sequence\t"
                                 "Count\t"
                                 "PValue\t"
                                 "Obs/Exp Max\t"
                                 "Max Obs/Exp Position\n";
  auto r = std::format(start_tag, grades.grade(label));
  r += header;
  std::ranges::for_each(
    results, [&](const auto &x) { r += std::format("{}\n", x.string()); });
  return r + end_module_tag;
}
