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

#include "falco_utils.hpp"
#include "falco_word.hpp"

#include <algorithm>
#include <cstdint>
#include <format>
#include <functional>
#include <limits>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

auto
duplication_results::operator+=(const duplication_results &rhs)
  -> const duplication_results & {
  for (const auto &[k, v] : rhs.dups)
    dups[k] += v;
  add_unique_seqs = add_unique_seqs && rhs.add_unique_seqs;
  n_unique = std::size(dups);
  limit_count += rhs.limit_count;
  return *this;
}

[[nodiscard]] auto
duplication_results::format_overrepresented(const std::uint64_t n_reads) const
  -> std::string {
  static constexpr auto start_tag = ">>Overrepresented sequences\t{}\n";
  static constexpr auto header = "#Sequence\t"
                                 "Count\t"
                                 "Percentage\t"
                                 "Possible Source\n";
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
      r += std::format("{}\t{}\t{}\n", seq, n_obs, pct(as_frac(n_obs, tot)));
  }
  return r + end_module_tag;
}

[[nodiscard]] auto
duplication_results::format_duplication_levels(
  [[maybe_unused]] const std::uint64_t n_reads) const -> std::string {
  static constexpr auto start_tag = ">>Sequence Duplication Levels\t{}\n"
                                    "#Total Deduplicated Percentage\t{:.6f}\n";
  static constexpr auto header = "#Duplication Level\t"
                                 "Percentage of deduplicated\t"
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

  const auto make_bins = [&](const auto &hist) {
    std::vector<std::uint64_t> binned(std::size(bin_breaks), 0);
    auto bin_id = 0ul;
    for (const auto [i, h] : std::views::enumerate(hist)) {
      // NOLINTNEXTLINE (*-pro-bounds-constant-array-index)
      if (bin_id < std::size(bin_breaks) && i >= bin_breaks[bin_id])
        ++bin_id;
      binned[bin_id] += h;
    }
    return binned;
  };

  [[maybe_unused]] const auto to_pct = [make_bins](const auto &v) {
    const auto sum = std::reduce(std::cbegin(v), std::cend(v));
    const auto p = [&](const auto d) { return pct(as_frac(d, sum)); };
    return make_bins(v) | std::views::transform(p) |
           std::ranges::to<std::vector>();
  };

  const auto binned_dedup = make_bins(hist_dedup);
  const auto binned_mass = make_bins(hist_mass);

  const auto grade = get_grade(grade_cutoffs, frac_dedup);
  auto r = std::format(start_tag, grade, pct(frac_dedup));
  r += header;
  auto not_first = 0;  // ADS: start output at 1 and not 0
  for (auto [label, dedup, mass] :
       std::views::zip(bin_labels, binned_dedup, binned_mass))
    if (not_first++)
      r += std::format("{}\t{}\t{}\n", label, dedup, mass);

  return r + end_module_tag;
}
