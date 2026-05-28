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

#include "duplication_results.hpp"

#include "contaminants.hpp"
#include "falco_utils.hpp"
#include "falco_word.hpp"

#include <algorithm>
#include <cstdint>
#include <format>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

std::int32_t duplication_results::read_skip = default_read_skip;

auto
duplication_results::initialize(const std::uint64_t est_n_reads) -> void {
  read_skip = est_n_reads < max_n_reads_total
                ? 0
                : static_cast<std::int32_t>(est_n_reads / max_n_reads_total);
}

auto
duplication_results::operator+=(const duplication_results &rhs)
  -> const duplication_results & {
  for (const auto &[k, v] : rhs.dups)
    dups[k] += v;
  return *this;
}

[[nodiscard]] auto
duplication_results::format_overrepresented() const -> std::string {
  static constexpr auto start_tag = ">>Overrepresented sequences\t{}\n";
  static constexpr auto header = "#Sequence\t"
                                 "Count\t"
                                 "Percentage\t"
                                 "Possible Source\n";
  const auto n_reads = std::reduce(std::cbegin(std::views::values(dups)),
                                   std::cend(std::views::values(dups)));
  const auto cutoff = static_cast<double>(n_reads) * overrep_cutoff;
  const auto max_n_obs = std::ranges::max(std::views::values(dups));
  const auto gt_cutoff = [&](const auto p) { return p.second >= cutoff; };
  const auto rev_p = [&](const auto p) { return std::pair{p.second, p.first}; };
  auto overrep = dups | std::views::filter(gt_cutoff) |
                 std::views::transform(rev_p) | std::ranges::to<std::vector>();
  std::ranges::sort(overrep, std::greater{});
  const auto reduce = [](const auto &v) {
    return std::reduce(std::cbegin(v), std::cend(v));
  };
  const auto tot = reduce(std::views::values(dups));
  const auto grade = get_grade(overrep_grade_cutoffs, as_frac(max_n_obs, tot));
  auto r = std::format(start_tag, grade);
  if (!overrep.empty()) {
    r += header;
    for (const auto &[n_obs, seq] : overrep)
      r += std::format("{}\t{}\t{}\t{}\n", seq, n_obs, pct(as_frac(n_obs, tot)),
                       match_contaminant(seq.string(), contaminants));
  }
  return r + end_module_tag;
}

[[nodiscard]] auto
duplication_results::format_duplication_levels() const -> std::string {
  static constexpr auto start_tag = ">>Sequence Duplication Levels\t{}\n"
                                    "#Total Deduplicated Percentage\t{:.6f}\n";
  static constexpr auto header = "#Duplication Level\t"
                                 "Percentage of total\n";
  // clang-format off
  static constexpr auto bin_breaks = std::array{
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    50,
    100,
    500,
    1'000,
    5'000,
    10'000,
    std::numeric_limits<int>::max(),
  };
  static constexpr auto bin_labels = std::array{
    "0",
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    ">10",
    ">50",
    ">100",
    ">500",
    ">1k",
    ">5k",
    ">10k",
  };
  static_assert(std::size(bin_labels) == std::size(bin_breaks));
  // clang-format on

  const auto max_dup = std::ranges::max(std::views::values(dups));

  std::vector<std::uint64_t> hist_dedup(max_dup + 1);
  for (const auto n_copies : std::views::values(dups))
    ++hist_dedup[n_copies];
  const auto dedup_sum =
    std::reduce(std::cbegin(hist_dedup), std::cend(hist_dedup));

  const auto hist_mass =
    std::views::transform(std::views::enumerate(hist_dedup), [](const auto x) {
      return std::get<0>(x) * std::get<1>(x);
    });
  const auto mass_sum =
    std::reduce(std::cbegin(hist_mass), std::cend(hist_mass));

  const auto frac_dedup = as_frac(dedup_sum, mass_sum);
  const auto grade = get_grade(grade_cutoffs, frac_dedup);
  auto r = std::format(start_tag, grade, pct(frac_dedup));
  r += header;

  const auto make_bins = [&](const auto &hist) {
    std::vector<std::uint64_t> binned(std::size(bin_breaks), 0);
    auto b_itr = std::cbegin(bin_breaks);
    for (const auto [i, h] : std::views::enumerate(hist)) {
      // ADS: clang-tidy false positive?
      // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-pointer-arithmetic)
      b_itr += (b_itr < std::cend(bin_breaks) && i >= *b_itr);
      binned[std::distance(std::cbegin(bin_breaks), b_itr)] += h;
    }
    return binned;
  };

  const auto to_pct = [make_bins](const auto &v) {
    const auto sum = std::reduce(std::cbegin(v), std::cend(v));
    const auto p = [&](const auto d) { return pct(as_frac(d, sum)); };
    return make_bins(v) | std::views::transform(p) |
           std::ranges::to<std::vector>();
  };

  const auto binned_mass = to_pct(hist_mass);
  for (const auto [label, mass] :
       std::views::zip(bin_labels, binned_mass) | std::views::drop(1))
    r += std::format("{}\t{:.6g}\n", label, mass);

  return r + end_module_tag;
}
