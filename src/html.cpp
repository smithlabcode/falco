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

#include "html.hpp"

#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "html_boilerplate.hpp"
#include "kmer_counter.hpp"  // IWYU pragma: keep

#define FMT_HEADER_ONLY
#include "fmt/base.h"
#include "fmt/chrono.h"  // IWYU pragma: keep
#include "fmt/format.h"
#include "fmt/ranges.h"

#include <config.h>

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <initializer_list>
#include <iterator>
#include <map>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

[[nodiscard]] auto
get_summary(const file_grades &grades) -> std::string {
  static constexpr auto summary_fmt = R"(
<div class="summary"><h2>Summary</h2>
<ul>
{}
</ul>
</div>
)";
  static constexpr auto item_fmt =
    R"(<li><a class="{}" href="#{}">{}</a></li>)";
  std::vector<std::string> sections;
  for (const auto &name : section_names)
    if (grades.is_configured(name))
      sections.emplace_back(fmt::format(item_fmt, grades.grade(name), name,
                                        grades.get_title(name)));
  return fmt::format(summary_fmt, fmt::join(sections, "\n"));
}

[[nodiscard]] auto
get_html_module(const std::string &label, const std::string &text,
                const file_grades &grades) -> std::string {
  static constexpr auto module_fmt =
    R"(<div class="module">
<h2 class="{}" id="{}">{}: {}</h2>
{}
</div>)";
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(module_fmt, label, grade, title, grade, text);
}

[[nodiscard]] auto
falco_get_html(const file_info &info, const file_grades &grades,
               const std::string &analysis_modules) -> std::string {
  return fmt::format(falco_html_body,
                     fmt::arg("date", std::chrono::system_clock::now()),
                     fmt::arg("filename", info.name),           //
                     fmt::arg("style", style),                  //
                     fmt::arg("summary", get_summary(grades)),  //
                     fmt::arg("modules", analysis_modules),     //
                     fmt::arg("version", VERSION));
}

[[nodiscard]] auto
sequence_length_html(const std::vector<std::uint64_t> &lengths,
                     const file_grades &grades) -> std::string {
  static constexpr auto label = "sequence_length";
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
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(html_module_fmt, grade, label, title, grade,
                     std::format(plot_fmt, lines_data));
}

[[nodiscard]] auto
gc_sequence_html(const falco::gc_content_array &gc_content,
                 const file_grades &grades) -> std::string {
  static constexpr auto label = "gc_sequence";
  static constexpr auto smoothing_window = 5;
  static constexpr auto plot_fmt =
    R"(<div id="gc_sequence_plot"></div>
<script>
Plotly.newPlot("gc_sequence_plot",
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
  const auto x = std::views::iota(1, std::ssize(gc_content) + 1);
  const auto total_count =
    std::reduce(std::cbegin(gc_content), std::cend(gc_content));
  const auto theor = get_theoretical_distribution(gc_content, total_count);
  const auto gc_smoothed = smooth_gc_content(gc_content, smoothing_window);
  const auto lines_data = fmt::format(gc_fmt, x, gc_smoothed, x, theor);

  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(html_module_fmt, grade, label, title, grade,
                     fmt::format(plot_fmt, lines_data));
}

[[nodiscard]] auto
sequence_html(const std::vector<falco::nuc_array> &nucs,
              const std::vector<base_group_t> &groups,
              const file_grades &grades) -> std::string {
  static constexpr auto label = "sequence";
  static constexpr auto plot_fmt =
    R"(<div id="sequence_plot"></div>
<script>
Plotly.newPlot("sequence_plot",
[{}],
{{
margin: {{t: 0}},
showlegend: true,
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
  static constexpr auto base_permutation = {0, 1, 3, 2};
  // ADS: the permutation is likely wrong...
  static constexpr auto bases = "ACTG";  // ACTG?
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
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(html_module_fmt, grade, label, title, grade,
                     fmt::format(plot_fmt, fmt::join(r, ",\n")));
}

[[nodiscard]] auto
n_content_html(const std::vector<std::uint64_t> &n_counts,
               const std::vector<falco::nuc_array> &nucs,
               const std::vector<base_group_t> &groups,
               const file_grades &grades) -> std::string {
  static constexpr auto label = "n_content";
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
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    fmt::format(
      plot_fmt, fmt::join(groups | std::views::transform(make_tag), ","),
      fmt::join(std::views::zip_transform(pct_for_pos, nucs, n_counts), ",")));
}

[[nodiscard]] auto
quality_sequence_html(const falco::qual_array &qual_by_read,
                      const file_grades &grades) -> std::string {
  static constexpr auto label = "quality_sequence";
  static constexpr auto plot_fmt =
    R"(<div id="quality_sequence_plot"></div>
<script>
Plotly.newPlot("quality_sequence_plot",
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

  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    fmt::format(plot_fmt, fmt::join(x, ","), fmt::join(y, ",")));
}

[[nodiscard]] auto
quality_base_html(const std::vector<falco::qual_array> &qual,
                  const std::vector<base_group_t> &groups,
                  const file_grades &grades) -> std::string {
  static constexpr auto label = "quality_base";
  static constexpr auto plot_fmt =
    R"(<div id="quality_base_plot"></div>
<script>
Plotly.newPlot("quality_base_plot",
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
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(html_module_fmt, grade, label, title, grade,
                     fmt::format(plot_fmt, fmt::join(lines, ",\n")));
}

[[nodiscard]] auto
basic_stats_html(const file_info &info, const std::uint64_t n_reads,
                 const std::uint64_t min_read_len,
                 const std::uint64_t max_read_len, const std::uint64_t total_gc,
                 const std::uint64_t total_nucs,
                 const file_grades &grades) -> std::string {
  static constexpr auto label = "basic_stats";
  static constexpr auto table_fmt =
    R"(<table><thead><tr><th>Measure</th><th>Value</th></tr></thead><tbody>
<tr><td>Filename</td><td>{filename_stem}</td></tr>
<tr><td>File type</td><td>{file_type}</td></tr>
<tr><td>Encoding</td><td>{encoding}</td></tr>
<tr><td>Total Sequences</td><td>{n_reads}</td></tr>
<tr><td>Sequences Flagged As Poor Quality</td><td>{n_poor}</td></tr>
<tr><td>Sequence length</td><td>{lengths_label}</td></tr>
<tr><td>%GC:</td><td>{gc_sequence_label:.1f}</td></tr>
</tbody></table>)";
  const auto lengths_label =
    min_read_len == max_read_len
      ? std::format("{}", max_read_len)
      : std::format("{}-{}", min_read_len, max_read_len);
  const auto gc_content_frac = pct(as_frac(total_gc, total_nucs));

  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    fmt::format(table_fmt, fmt::arg("filename_stem", info.name),
                fmt::arg("file_type", info.description),
                fmt::arg("encoding", to_string(info.encoding)),
                fmt::arg("n_reads", n_reads),
                fmt::arg("n_poor", 0),  // ADS: where to calcualte this?
                fmt::arg("lengths_label", lengths_label),
                fmt::arg("gc_sequence_label", gc_content_frac)));
}

[[nodiscard]] auto
tile_html(const tile_processor::tiles_centered_t &centered,
          const std::vector<base_group_t> &groups,
          const file_grades &grades) -> std::string {
  static constexpr auto label = "tile";
  static constexpr auto n_quants = 20.0;
  // ADS: ??? (-10: red, 0: light blue, +10: dark blue)
  static constexpr auto tiles_plot_fmt =
    R"""([{{
x: [{}],
y: [{}],
z: [{}],
type: "heatmap",
colorscale: [
[0.0, "rgb(210,65,83)"],
[{}, "rgb(178,236,254)"],
[1.0, "rgb(34,57,212)"]
],
showscale: true,
}}],
{{
margin: {{t: 0}},
showlegend: false,
xaxis: {{title: "Base position"}},
yaxis: {{title: "tile", type: "category"}},
}}
)""";
  static constexpr auto plot_fmt = R"(<div id="{}"></div>
<script>Plotly.newPlot("{}",
{}
);</script>
)";
  assert(get_max_size(centered) == std::size(groups) || centered.empty());
  if (centered.empty())
    return {};  // ADS: in case we are here but tile analysis was not done
  const auto tag = [&](const auto &g) { return make_group_tag_quoted(g); };

  const auto [minval, maxval] =
    std::ranges::minmax(centered | std::views::values | std::views::join);
  // discretize quantiles so plotly understands color scheme (???)
  const auto midpoint = minval / (minval - maxval);
  const auto mid_discrete = std::round(n_quants * midpoint) / n_quants;

  // z: one formatted vector of quality z scores per tile
  const auto format1 = [](const auto &c) {
    // ADS: note the format and precision (3) here
    return fmt::format("[{:.3f}]", fmt::join(c, ","));  // inner lists
  };
  const auto z = std::views::transform(centered | std::views::values, format1);

  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    fmt::format(
      plot_fmt, "tiles_plot", "tiles_plot",
      fmt::format(tiles_plot_fmt,
                  fmt::join(groups | std::views::transform(tag), ","),  // x
                  fmt::join(centered | std::views::keys, ","),          // y
                  fmt::join(z, ","),                                    // z
                  mid_discrete)));
}

[[nodiscard]] auto
kmer_html(const std::vector<kmer_result> &results,
          const file_grades &grades) -> std::string {
  static constexpr auto label = "kmer";
  static constexpr auto plot_format = R"(<div id="{}"></div>
<script>Plotly.newPlot("{}",
{}
);</script>
)";
  static constexpr auto kmer_line_format =
    R"({{
x: [{}],
y: [{}],
type: "line",
name: "{}",
}})";
  const auto get_pos = [&](const auto &x) { return x.pos; };
  const auto xlim = std::ranges::max(std::views::transform(results, get_pos));
  auto mostly0 = std::vector(xlim, 0.0);  // change one position each kmer
  const auto xvals = std::views::iota(0u, xlim);
  auto prev = 0;
  const auto format1 = [&](const auto &k) {
    // x: positions in read; y: zeros except at 'pos' for the kmer
    mostly0[prev] = 0.0;  // replace the zero for previous position 'k.pos'
    mostly0[k.pos] = std::log2(k.obs_exp);
    prev = k.pos;
    return fmt::format(kmer_line_format, fmt::join(xvals, ", "),
                       fmt::join(mostly0, ", "), k.decode());
  };
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    fmt::format(
      plot_format, "kmer_plot", "kmer_plot",
      fmt::format("[{}]",
                  fmt::join(std::views::transform(results, format1), ",\n"))));
}
