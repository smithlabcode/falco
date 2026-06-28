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
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "falco_word.hpp"
#include "html.hpp"
#include "run_mode.hpp"

#include "boost/boost_unordered.hpp"

#define FMT_HEADER_ONLY
#include "fmt/format.h"
#include "fmt/ranges.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <format>
#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// ADS: Number of bins is not user-adjustable.

static constexpr auto n_bins = 16;

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

// clang-format off

// ADS: previously for plots: std::array{"1", "2", ..., "5k+", "10k+"}
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

// clang-format on

// ADS: One extra bin (at max) for analyis, and one more (at min) among labels
static_assert(std::size(bin_labels) == std::size(bin_breaks) &&
              n_bins + 1 == std::size(bin_labels));

auto
duplication_results::get_n_reads() const -> std::uint64_t {
  const auto dups_v = std::views::values(dups);
  return std::reduce(std::cbegin(dups_v), std::cend(dups_v));
}

[[nodiscard]] auto
duplication_results::get_overrepresented(const std::uint64_t n_reads) const
  -> std::vector<overrep_t> {
  const auto cutoff = static_cast<double>(n_reads) * overrep_cutoff;
  const auto gt_cutoff = [&](const auto p) { return p.second >= cutoff; };
  const auto rev_p = [&](const auto p) { return std::pair{p.second, p.first}; };
  auto overrep = dups | std::views::filter(gt_cutoff) |
                 std::views::transform(rev_p) | std::ranges::to<std::vector>();
  std::ranges::sort(overrep, std::greater{});
  std::vector<overrep_t> ret;
  for (const auto &[n_obs, seq] : overrep)
    ret.emplace_back(seq, n_obs, pct(as_frac(n_obs, n_reads)),
                     match_contaminant(seq.string()));
  return ret;
}

auto
duplication_results::initialize(const run_mode &mode,
                                const file_info &info) -> void {
  read_skip =
    info.n_reads_est < max_n_reads_total
      ? 0
      : static_cast<std::int32_t>(info.n_reads_est / max_n_reads_total);
  if (!mode.do_dups()) {
    // ADS: disabling dups analysis
#ifdef ORIGINAL_DUPS
    max_reads_to_hash = 0;
#else  // NOT ORIGINAL_DUPS
    read_idx = std::numeric_limits<std::int64_t>::max();
#endif
  }
}

auto
duplication_results::operator+=(const duplication_results &rhs)
  -> const duplication_results & {
  for (const auto &[k, v] : rhs.dups)
    dups[k] += v;
  return *this;
}

[[nodiscard]] auto
get_grade_overrepresented(const std::uint64_t n_reads,
                          const duplication_results &dr) -> std::string {
  static constexpr auto label = "overrepresented";
  const auto max_n_obs =
    dr.dups.empty() ? 0LU : std::ranges::max(std::views::values(dr.dups));
  return grader_set::get_grade(label, as_frac(max_n_obs, n_reads));
}

[[nodiscard]] auto
overrepresented_report(const std::vector<overrep_t> &overrep,
                       const file_grades &grades) -> std::string {
  static constexpr auto label = "overrepresented";
  static constexpr auto start_tag = ">>Overrepresented sequences\t{}\n";
  static constexpr auto header = "#Sequence\t"
                                 "Count\t"
                                 "Percentage\t"
                                 "Possible Source\n";
  auto r = std::format(start_tag, grades.grade(label));
  if (!overrep.empty()) {
    r += header;
    for (const auto &[seq, n_obs, pct_val, contam_id] : overrep)
      r += std::format("{}\t{}\t{:.3g}\t{}\n", seq, n_obs, pct_val,
                       get_contam_name(contam_id));
  }
  return r + end_module_tag;
}

[[nodiscard]] static auto
make_bins(const auto &breaks, const auto &hist) {
  std::vector<std::uint64_t> binned(std::size(breaks), 0);
  auto b_itr = std::cbegin(breaks);
  for (const auto [i, h] : std::views::enumerate(hist)) {
    // ADS: clang-tidy false positive?
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    b_itr += (b_itr < std::cend(breaks) && i >= *b_itr);
    binned[std::distance(std::cbegin(breaks), b_itr)] += h;
  }
  return binned;
}

[[nodiscard]] auto
duplication_results::get_dups_summary() const -> dup_summary_t {
  if (dups.empty())
    return {};
  const auto max_dup = std::ranges::max(std::views::values(dups));
  std::vector<std::uint64_t> hist_dedup(max_dup + 1);
  for (const auto n_copies : std::views::values(dups))
    ++hist_dedup[n_copies];
  const auto hist_mass =
    std::views::transform(
      std::views::enumerate(hist_dedup),
      [](const auto x) { return std::get<0>(x) * std::get<1>(x); }) |
    std::ranges::to<std::vector>();
  // ADS: move needed to avoid copying during creation of the tuple
  return dup_summary_t{
    max_dup,
    get_n_reads(),
    std::move(hist_mass),
    std::move(hist_dedup),
  };
}

[[nodiscard]] auto
get_grade_duplication(const dup_summary_t &summary) -> std::string {
  static constexpr auto label = "duplication";
  const auto reduce = [](const auto &v) {
    return std::reduce(std::cbegin(v), std::cend(v));
  };
  const auto frac_dedup =
    as_frac(reduce(summary.hist_dedup), reduce(summary.hist_mass));
  return grader_set::get_grade(label, frac_dedup);
}

[[nodiscard]] auto
duplication_report(const dup_summary_t &summary,
                   const file_grades &grades) -> std::string {
  static constexpr auto label = "duplication";
  static constexpr auto start_tag = ">>Sequence Duplication Levels\t{}\n"
                                    "#Total Deduplicated Percentage\t{:.6f}\n";
  static constexpr auto header = "#Duplication Level\t"
                                 "Percentage of total\n";
  const auto to_pct = [&](const auto &v) {
    const auto sum = std::reduce(std::cbegin(v), std::cend(v));
    const auto p = [&](const auto d) { return pct(as_frac(d, sum)); };
    return make_bins(bin_breaks, v) | std::views::transform(p) |
           std::ranges::to<std::vector>();
  };
  const auto binned_mass_pct = to_pct(summary.hist_mass);
  const auto reduce = [](const auto &v) {
    return std::reduce(std::cbegin(v), std::cend(v));
  };
  const auto frac_dedup =
    as_frac(reduce(summary.hist_dedup), reduce(summary.hist_mass));
  auto r = std::format(start_tag, grades.grade(label), pct(frac_dedup));
  r += header;
  for (const auto [label, mass] :
       std::views::zip(bin_labels, binned_mass_pct) | std::views::drop(1))
    r += std::format("{}\t{:.3g}\n", label, mass);  // Percentage format is .3g
  return r + end_module_tag;
}

[[nodiscard]] auto
overrepresented_html(const std::vector<overrep_t> &overrep,
                     const file_grades &grades) -> std::string {
  static constexpr auto label = "overrepresented";
  static constexpr auto html_table = R"(<table>
<thead>
<tr>
<th>Sequence</th>
<th>Count</th>
<th>Percentage</th>
<th>Possible Source</th>
</tr>
</thead>
<tbody>
{}
</tbody>
</table>
)";
  static constexpr auto html_table_row_fmt =  // Percentage format is .3g
    "<tr><td>{}</td><td>{}</td><td>{:.3g}</td><td>{}</td></tr>";
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  if (overrep.empty())
    return fmt::format(html_module_fmt, grade, label, title, grade,
                       "No overrepresented sequences");
  std::vector<std::string> rows;
  for (const auto &[seq, n_obs, pct_val, contam_id] : overrep)
    rows.push_back(fmt::format(html_table_row_fmt, seq.string(), n_obs, pct_val,
                               get_contam_name(contam_id)));
  return fmt::format(html_module_fmt, grade, label, title, grade,
                     fmt::format(html_table, fmt::join(rows, "\n")));
}

[[nodiscard]] auto
duplication_html(const dup_summary_t &summary,
                 const file_grades &grades) -> std::string {
  static constexpr auto label = "duplication";
  static constexpr auto plot_format = R"(<div id="duplication_plot"></div>
<script>Plotly.newPlot("duplication_plot",
[{{
x: {},
y: {},
type: "line",
line: {{color: "blue"}},
name: "total sequences",
}},
{{
x: {},
y: {},
type: "line",
line: {{color: "red"}},
name: "deduplicated sequences",
}}
],
{{
margin: {{t: 0}},
showlegend: true,
xaxis: {{
title: "Duplication rate",
tickvals: {},
ticktext: {},
}},
yaxis: {{title: "% of sequences"}},
}}
);</script>
)";
  static constexpr auto x = std::views::iota(1, n_bins + 1);
  static constexpr auto x_text =  // has one extra element at "0"
    std::ranges::subrange(std::cbegin(bin_labels), std::cend(bin_labels));
  const auto to_pct = [&](const auto &v) {
    const auto sum = std::reduce(std::cbegin(v), std::cend(v));
    const auto p = [&](const auto d) { return pct(as_frac(d, sum)); };
    return make_bins(bin_breaks, v) | std::views::transform(p) |
           std::ranges::to<std::vector>();
  };
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    fmt::format(plot_format,                                          //
                x, to_pct(summary.hist_mass) | std::views::drop(1),   // y_tot,
                x, to_pct(summary.hist_dedup) | std::views::drop(1),  // y_dedup
                x, x_text | std::views::drop(1)  // tickvals
                ));
}
