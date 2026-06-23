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

#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "quality_score.hpp"

#include <config.h>

#define FMT_HEADER_ONLY
#include "fmt/base.h"
#include "fmt/format.h"
#include "fmt/ranges.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <format>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <limits>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

[[nodiscard]] auto
get_grade_read_lengths(const std::vector<std::uint64_t> &lengths)
  -> std::string {
  using std::literals::string_literals::operator""s;
  const bool has_empty_reads = std::size(lengths) > 0 && lengths[0] > 0;
  const auto n_lengths =
    std::ranges::count_if(lengths, [](const auto x) { return x > 0; });
  return has_empty_reads ? "fail"s : (n_lengths > 1 ? "warn"s : "pass"s);
}

[[nodiscard]] auto
format_read_lengths(const std::vector<std::uint64_t> &lengths,
                    const std::string &grade) -> std::string {
  static constexpr auto start_tag = ">>Sequence Length Distribution\t{}\n";
  static constexpr auto header = "#Length\t"
                                 "Count\n";
  const auto eq0 = [](const auto &x) { return std::get<1>(x) == 0; };
  auto r = std::format(start_tag, grade);
  r += header;
  auto to_report = std::views::enumerate(lengths) | std::views::drop_while(eq0);
  for (const auto [idx, len] : to_report)
    r += std::format("{}\t{}\n", idx, len);
  r += end_module_tag;
  return r;
}

[[nodiscard]] static inline auto
get_theoretical_distribution(const auto &gc, const auto &total_count) {
  static constexpr auto mode_width = 0.90;  // ADS: where does this come from?
  // calculate deviation of a hist from a "normal" with same mode and sd
  const auto n_bins = std::size(gc);
  const auto gc_beg = std::cbegin(gc);
  if (total_count <= 1)  // we will be dividing by (total_count - 1)
    return std::vector<double>(std::size(gc));

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
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
    return std::exp(-cntr_sq(val) / (2.0 * sd * sd));
  };
  auto normed = std::views::iota(0u, n_bins) |
                std::views::transform(to_normal) |
                std::ranges::to<std::vector>();
  const auto denom = std::reduce(std::cbegin(normed), std::cend(normed));
  std::ranges::transform(normed, std::begin(normed),
                         [&](const auto x) { return x * total_count / denom; });
  return normed;
}

[[nodiscard]] static inline auto
smooth(const auto &data, const auto window_size) {
  const auto get_mean = [&](std::ranges::viewable_range auto &&r) {
    return as_frac(std::reduce(std::cbegin(r), std::cend(r)), std::size(r));
  };
  assert(window_size < std::size(data));
  std::vector<double> smoothed;
  for (auto w = 1; w < (window_size + 1) / 2; ++w)
    smoothed.push_back(get_mean(
      std::ranges::subrange(std::cbegin(data), std::cbegin(data) + w)));
  for (const auto &window : data | std::views::slide(window_size))
    smoothed.push_back(get_mean(window));
  for (auto w = (window_size + 1) / 2; w > 1; --w)
    smoothed.push_back(get_mean(
      std::ranges::subrange(std::cend(data) - w + 1, std::cend(data))));
  return smoothed;
}

[[nodiscard]] static inline auto
sum_deviation_from_normal(const auto &gc) {
  const auto gc_beg = std::cbegin(gc);
  const auto gc_end = std::cend(gc);
  const auto total_count = std::reduce(gc_beg, gc_end);
  if (total_count <= 1)  // we will be dividing by (total_count - 1)
    return 0.0;

  const auto theor = get_theoretical_distribution(gc, total_count);
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
  return identify_grade(grade_cutoffs, deviation);
}

[[nodiscard]] auto
format_gc_content(const falco::gc_content_array &gc_content,
                  const std::string &grade) -> std::string {
  static constexpr auto start_tag = ">>Per sequence GC content\t{}\n";
  static constexpr auto header = "#GC Content\t"
                                 "Count\n";
  auto r = std::format(start_tag, grade);
  r += header;
  for (const auto [idx, gc] : std::views::enumerate(gc_content))
    r += std::format("{}\t{}\n", idx, gc);
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
get_grade_base_comp(const std::vector<falco::nuc_array> &nucs) -> std::string {
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
  return identify_grade(grade_cutoffs, max_diff);
}

[[nodiscard]] auto
format_base_comp(const std::vector<falco::nuc_array> &nucs,
                 const std::vector<base_group_t> &groups,
                 const std::string &grade) -> std::string {
  static constexpr auto start_tag = ">>Per base sequence content\t{}\n";
  static constexpr auto header = "#Base\t"
                                 "G\tA\tT\tC\n";
  static constexpr auto base_permutation = {3, 0, 2, 1};
  auto r = std::format(start_tag, grade);
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
get_grade_n_content(const std::vector<std::uint64_t> &n_counts,
                    const std::vector<falco::nuc_array> &nucs) -> std::string {
  // ADS: at this point 'nucs' should not include counts of 'N'
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.05, "pass"},
    std::pair{0.20, "warn"},
    std::pair{1.00, "fail"},
  };
  const auto max_idx =
    std::distance(std::cbegin(n_counts), std::ranges::max_element(n_counts));
  const auto total_non_n =
    std::reduce(std::cbegin(nucs[max_idx]), std::cend(nucs[max_idx]));
  const auto max_n_frac =
    as_frac(n_counts[max_idx], n_counts[max_idx] + total_non_n);
  return identify_grade(grade_cutoffs, max_n_frac);
}

[[nodiscard]] auto
format_n_content(const std::vector<std::uint64_t> &n_counts,
                 const std::vector<falco::nuc_array> &nucs,
                 const std::vector<base_group_t> &groups,
                 const std::string &grade) -> std::string {
  static constexpr auto start_tag = ">>Per base N content\t{}\n";
  static constexpr auto header = "#Base\t"
                                 "N-Count\n";
  auto r = std::format(start_tag, grade);
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
  return identify_grade(grade_cutoffs, qual_val_mode);
}

[[nodiscard]] auto
format_qual_by_read(const falco::qual_array &qual_by_read,
                    const std::string &grade) -> std::string {
  static constexpr auto start_tag = ">>Per sequence quality scores\t{}\n";
  static constexpr auto header = "#Quality\t"
                                 "Count\n";
  auto r = std::format(start_tag, grade);
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
get_grade_qual_by_pos(const std::vector<falco::qual_array> &qual)
  -> std::string {
  // NOLINTBEGIN(cppcoreguidelines-avoid-magic-numbers)
  static constexpr auto lower_quartile = [](const auto &q) { return q[1]; };
  static constexpr auto get_median = [](const auto &q) { return q[2]; };
  // NOLINTEND(cppcoreguidelines-avoid-magic-numbers)
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
  return identify_grade(grade_cutoffs_median, min_qual_median);
}

[[nodiscard]] auto
format_qual_by_pos(const std::vector<falco::qual_array> &qual,
                   const std::vector<base_group_t> &groups,
                   const std::string &grade) -> std::string {
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
  auto r = std::format(start_tag, grade);
  r += header;
  for (const auto [idx, q] : std::views::enumerate(qual))
    r += std::format("{}\t{:.{}g}{}\n", make_group_tag(groups[idx]),
                     mean_tabular(q), digits, tab_sep(five_quants(q)));
  r += end_module_tag;
  return r;
}

[[nodiscard]] auto
get_grade_basic_stats() -> std::string {
  static constexpr auto default_grade = "pass";  // always a pass
  return default_grade;
}

[[nodiscard]] auto
format_basic_stats(const file_info &info, const std::uint64_t n_reads,
                   const std::uint64_t min_read_len,
                   const std::uint64_t max_read_len,
                   const std::uint64_t total_gc, const std::uint64_t total_nucs,
                   const std::string &grade) -> std::string {
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
  auto r = std::format(start_tag, VERSION, grade);
  r += header;
  r += std::format("Filename\t{}\n", info.name);
  r += std::format("File type\t{}\n", info.description);
  r += std::format("Encoding\t{}\n", to_string(info.encoding));
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

[[nodiscard]] auto
format_read_lengths_html(const std::vector<std::uint64_t> &lengths)
  -> std::string {
  static constexpr auto plot_fmt =
    R"(<div id="length_plot"></div>
<script>
Plotly.newPlot("length_plot",
[{}],
{{
margin: {{t: 0}},
showlegend: true,
xaxis: {{title: "Sequence length"}},
yaxis: {{title: "Number of sequences"}},
}});
</script>
)";
  static constexpr auto lengths_fmt = R"""({{
x: [{}],
y: [{}],
text: [{}],
type: "bar",
marker: {{color: "rgba(55,128,191,1.0)"}},
line: {{width : 2}},
name: "Sequence length distribution"
}})""";
  const auto length_values = std::views::iota(0LU, std::size(lengths));
  const auto add_bp = [&](const auto l) {
    return fmt::format(R"("{} bp")", l);
  };
  const auto lines_data = fmt::format(
    lengths_fmt, fmt::join(std::views::transform(length_values, add_bp), ","),
    fmt::join(lengths, ","), fmt::join(length_values, ","));
  return std::format(plot_fmt, lines_data);
}

[[nodiscard]] auto
format_gc_content_html(const falco::gc_content_array &gc_content)
  -> std::string {
  static constexpr auto smoothing_window = 5;
  static constexpr auto plot_fmt =
    R"(<div id="gc_content_plot"></div>
<script>
Plotly.newPlot("gc_content_plot",
[{}],
{{
margin: {{t: 0}},
showlegend: true,
xaxis: {{title: "% GC"}},
yaxis: {{title: "Density"}}
}});
</script>
)";
  static constexpr auto gc_fmt = R"({{
x: {},
y: {},
type: "line",
line: {{color: "red"}},
name: "GC distribution"
}},
{{
x: {},
y: {},
type: "line",
line: {{color : "blue"}},
name: "Theoretical distribution"
}})";
  const auto x = std::views::iota(0LU, std::size(gc_content));
  const auto total_count =
    std::reduce(std::cbegin(gc_content), std::cend(gc_content));
  const auto theor = get_theoretical_distribution(gc_content, total_count);
  const auto gc_smoothed = smooth(gc_content, smoothing_window);
  const auto lines_data = fmt::format(gc_fmt, x, gc_smoothed, x, theor);
  return fmt::format(plot_fmt, lines_data);
}

[[nodiscard]] auto
format_base_comp_html(const std::vector<falco::nuc_array> &nucs,
                      const std::vector<base_group_t> &groups) -> std::string {
  static constexpr auto plot_fmt =
    R"(<div id="base_comp_plot"></div>
<script>
Plotly.newPlot("base_comp_plot",
[{}],
{{
margin: {{t: 0}},
showlegend: false,
xaxis: {{title: "Base position"}},
yaxis: {{title: "Phread quality"}},
}});
</script>
)";
  static constexpr auto per_base_fmt = R"({{
x: [{}],
y: [{:.3f}],
mode: "lines",
name: "{}",
line: {{color : "{}"}}
}})";
  static constexpr auto base_permutation = {3, 0, 2, 1};
  // ADS: the permutation is likely wrong...
  static constexpr auto bases = "ATGC";  // ATGC
  static constexpr auto base_to_color = std::array{
    "green",
    "blue",
    "red",
    "black",
  };
  const auto sum = [&](const auto &nucs_by_pos) {
    return std::reduce(std::cbegin(nucs_by_pos), std::cend(nucs_by_pos));
  };
  assert(std::size(nucs) == std::size(groups));
  const auto x = groups | std::views::transform([&](const auto &g) {
                   return make_group_tag_quoted(g);
                 });
  const auto total_by_pos = nucs | std::views::transform(sum);
  std::vector<std::string> r;
  // NOLINTBEGIN(*-constant-array-index,*-pointer-arithmetic)
  for (const auto idx : base_permutation) {
    const auto pct_for_pos = [idx](const auto &nucs_for_pos, const auto tot) {
      return pct(as_frac(nucs_for_pos[idx], tot));
    };
    const auto y = std::views::zip_transform(pct_for_pos, nucs, total_by_pos);
    r.emplace_back(fmt::format(per_base_fmt, fmt::join(x, ","),
                               fmt::join(y, ","), bases[idx],
                               base_to_color[idx]));
  }
  // NOLINTEND(*-constant-array-index,*-pointer-arithmetic)
  return fmt::format(plot_fmt, fmt::join(r, ",\n"));
}

[[nodiscard]] auto
format_n_content_html(const std::vector<std::uint64_t> &n_counts,
                      const std::vector<falco::nuc_array> &nucs,
                      const std::vector<base_group_t> &groups) -> std::string {
  static constexpr auto plot_fmt =
    R"(<div id="n_content_plot"></div>
<script>
Plotly.newPlot("n_content_plot",
[{{
x: [{}],
y: [{:.6g}],
type: "line",
line: {{color : "red"}},
name: "Fraction of N reads per base"
}}],
{{
margin: {{t: 0}},
showlegend: true,
xaxis: {{title: "Base position"}},
yaxis: {{title: "% N"}},
}}
);</script>
)";
  const auto pct_for_pos = [](const auto &nucs_by_pos, const auto n_count) {
    return pct(as_frac(
      n_count, std::reduce(std::cbegin(nucs_by_pos), std::cend(nucs_by_pos))));
  };
  const auto make_tag = [&](const auto &g) { return make_group_tag_quoted(g); };
  assert(std::size(nucs) == std::size(groups));
  return fmt::format(
    plot_fmt, fmt::join(groups | std::views::transform(make_tag), ","),
    fmt::join(std::views::zip_transform(pct_for_pos, nucs, n_counts), ","));
}

[[nodiscard]] auto
format_qual_by_read_html(const falco::qual_array &qual_by_read) -> std::string {
  static constexpr auto plot_fmt =
    R"(<div id="qual_by_read_plot"></div>
<script>
Plotly.newPlot("qual_by_read_plot",
[{{
x: [{}],
y: [{}],
type: "line",
line: {{color: "red"}},
name: "Sequence quality distribution"
}}],
{{
margin: {{t: 0}},
showlegend: true,
xaxis: {{title: "Phread quality"}},
yaxis: {{title: "Density"}},
}});
</script>
)";
  // output quality values between first non-zero and last zero
  const auto gt0 = [&](const auto x) { return x > 0; };
  const auto first_obs_itr = std::ranges::find_if(qual_by_read, gt0);
  const auto q_beg = std::cbegin(qual_by_read);
  const std::int64_t first_obs = std::distance(q_beg, first_obs_itr);
  const auto last_obs_subrange = std::ranges::find_last_if(qual_by_read, gt0);
  const std::int64_t last_obs =
    std::ssize(qual_by_read) - std::ssize(last_obs_subrange) + 1;
  assert(first_obs >= 0 && last_obs <= falco::max_qual_val);

  const auto x = std::views::iota(first_obs, last_obs);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const auto y = std::ranges::subrange{q_beg + first_obs, q_beg + last_obs};
  return fmt::format(plot_fmt, fmt::join(x, ","), fmt::join(y, ","));
}

[[nodiscard]] auto
format_qual_by_pos_html(const std::vector<falco::qual_array> &qual,
                        const std::vector<base_group_t> &groups)
  -> std::string {
  static constexpr auto plot_fmt =
    R"(<div id="qual_by_pos_plot"></div>
<script>
Plotly.newPlot("qual_by_pos_plot",
[{}],
{{
margin: {{t: 0}},
showlegend: false,
xaxis: {{title: "Base position"}},
yaxis: {{title: "Phread quality"}},
}});
</script>
)";
  static constexpr auto row_fmt =
    R"({{y: [{}], type: "box", name: "{}bp", marker: {{color : "{}"}}}})";
  const auto get_color = [&](const auto &q) {
    static constexpr auto median_error = 0.20;
    static constexpr auto median_warn = 0.25;
    static constexpr auto lquartile_error = 0.05;
    static constexpr auto lquartile_warn = 0.10;
    if (median_val(q) < median_error || lquart_val(q) < lquartile_error)
      return "red";
    if (median_val(q) < median_warn || lquart_val(q) < lquartile_warn)
      return "yellow";
    return "green";
  };
  std::vector<std::string> lines;
  for (const auto [idx, q] : std::views::enumerate(qual)) {
    const auto fq = five_quants(q);
    lines.emplace_back(fmt::format(row_fmt, fmt::join(fq, ", "),
                                   make_group_tag(groups[idx]), get_color(fq)));
  }
  return fmt::format(plot_fmt, fmt::join(lines, ",\n"));
}

[[nodiscard]] auto
format_basic_stats_html(const file_info &info, const std::uint64_t n_reads,
                        const std::uint64_t min_read_len,
                        const std::uint64_t max_read_len,
                        const std::uint64_t total_gc,
                        const std::uint64_t total_nucs) -> std::string {
  static constexpr auto table_fmt =
    R"(<table><thead><tr><th>Measure</th><th>Value</th></tr></thead><tbody>
<tr><td>Filename</td><td>{filename_stem}</td></tr>
<tr><td>File type</td><td>{file_type}</td></tr>
<tr><td>Encoding</td><td>{encoding}</td></tr>
<tr><td>Total Sequences</td><td>{n_reads}</td></tr>
<tr><td>Sequences Flagged As Poor Quality</td><td>{n_poor}</td></tr>
<tr><td>Sequence length</td><td>{lengths_label}</td></tr>
<tr><td>%GC:</td><td>{gc_content_label:.1f}</td></tr>
</tbody></table>)";
  const auto lengths_label =
    min_read_len == max_read_len
      ? std::format("{}", max_read_len)
      : std::format("{}-{}", min_read_len, max_read_len);
  const auto gc_content_frac = pct(as_frac(total_gc, total_nucs));
  return fmt::format(table_fmt, fmt::arg("filename_stem", info.name),
                     fmt::arg("file_type", info.description),
                     fmt::arg("encoding", to_string(info.encoding)),
                     fmt::arg("n_reads", n_reads),
                     fmt::arg("n_poor", 0),  // ADS: where to calcualte this?
                     fmt::arg("lengths_label", lengths_label),
                     fmt::arg("gc_content_label", gc_content_frac));
}
