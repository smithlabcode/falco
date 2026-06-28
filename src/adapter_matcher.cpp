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

#include "adapter_matcher.hpp"

#include "adapter_set.hpp"
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "html.hpp"
#include "run_mode.hpp"

#define FMT_HEADER_ONLY
#include "fmt/format.h"
#include "fmt/ranges.h"  // IWYU pragma: keep

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <format>
#include <functional>
#include <iterator>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <vector>

adapter_matcher::adapter_matcher() {
  static constexpr auto shift_plus = [&](const auto r, const auto c) {
    return (r << nibble_size) + encode_nibble(c);
  };
  static constexpr auto do_encoding = [&](const auto &a) {
    return std::accumulate(std::cbegin(a), std::cend(a), 0ul, shift_plus);
  };
  const auto &as = adapter_set::instance();
  assert(!as.instance().adapters.empty());
  n_adapters = adapter_set::n_adapters();
  adapter_size = adapter_set::adapter_size();
  std::ranges::transform(as.adapters, std::back_inserter(encoded_adapters),
                         do_encoding);
}

auto
adapter_matcher::operator+=(const adapter_matcher &rhs)
  -> const adapter_matcher & {
  two_dim_add(adap_counts, rhs.adap_counts);
  return *this;
}

[[nodiscard]] auto
adapter_matcher::get_grade(const std::uint64_t n_reads) const -> std::string {
  static constexpr auto label = "adapter";
  const auto max_count =
    std::ranges::max(adap_counts | std::views::transform(std::ranges::max));
  return grader_set::get_grade(label, as_frac(max_count, n_reads));
}

[[nodiscard]] auto
adapter_matcher::report(const std::uint64_t n_reads,
                        const std::vector<base_group_t> &groups,
                        const file_grades &grades) const -> std::string {
  static constexpr auto label = "adapter";
  static constexpr auto start_module_tag = ">>Adapter Content\t{}\n";
  static constexpr auto header = "#Position";
  static constexpr auto to_flat = [](const auto &data, const auto &fun) {
    return data | std::views::transform(fun) | std::views::join |
           std::ranges::to<std::string>();
  };
  auto r = std::format(start_module_tag, grades.grade(label));
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  const auto &adapter_names = adapter_set::instance().adapter_names;
  r += header + to_flat(adapter_names,
                        [](const auto &x) { return std::format("\t{}", x); });
  r += '\n';

  auto cumulative = adap_counts;
  for (auto [prev, curr] : cumulative | std::views::pairwise)
    std::ranges::transform(curr, prev, std::begin(curr), std::plus{});

  const auto n_pos = std::size(adap_counts);
  const auto lim = n_pos + 1 >= adapter_size ? n_pos - adapter_size + 1 : n_pos;
  cumulative.resize(lim);
  const auto fmt_pct_of_reads = [n_reads](const auto c) {
    return std::format("\t{:.6g}", pct(as_frac(c, n_reads)));
  };
  for (const auto [idx, cumul] : std::views::enumerate(cumulative))
    r += std::format("{}{}\n", make_group_tag(groups[idx]),
                     to_flat(cumul, fmt_pct_of_reads));
  return r + end_module_tag;
}

[[nodiscard]] auto
adapter_matcher::html(const std::uint64_t n_reads,
                      const std::vector<base_group_t> &groups,
                      const file_grades &grades) const -> std::string {
  static constexpr auto label = "adapter";
  static constexpr auto plot_fmt =
    R"(<div id="adapters_plot"></div>
<script>
Plotly.newPlot("adapters_plot",
{},
{{
margin: {{t: 0}},
showlegend: true,
xaxis: {{title: "Base position"}},
yaxis: {{title: "% sequences with adapter before position"}},
}});
</script>)";
  static constexpr auto adapter_fmt = R"({{
x: {},
y: {},
type: "line",
name: "{}",
}})";
  assert(std::size(groups) == std::size(adap_counts));
  // calcualte the x axis first
  const auto x = groups | std::views::transform([&](const auto &g) {
                   return make_group_tag_quoted(g);
                 });
  auto cumulative = adap_counts;
  for (auto [prev, curr] : cumulative | std::views::pairwise)
    std::ranges::transform(curr, prev, std::begin(curr), std::plus{});

  const auto n_pos = std::size(adap_counts);
  const auto lim = n_pos + 1 >= adapter_size ? n_pos - adapter_size + 1 : n_pos;
  cumulative.resize(lim);
  const auto pct_of_reads = [n_reads](const auto c) {
    return pct(as_frac(c, n_reads));
  };
  std::vector<std::string> html_by_adapter;
  const auto adapter_names = adapter_set::instance().adapter_names;
  for (const auto [adap_id, adap_name] : std::views::enumerate(adapter_names)) {
    const auto make_y = [&](const auto &cumul) {
      return pct_of_reads(cumul[adap_id]);
    };
    const auto y = cumulative | std::views::transform(make_y);
    html_by_adapter.emplace_back(fmt::format(adapter_fmt, x, y, adap_name));
  }
  const auto grade = grades.grade(label);
  const auto title = grades.get_title(label);
  return fmt::format(
    html_module_fmt, grade, label, title, grade,
    std::format(plot_fmt, fmt::format("[{:n:}]", html_by_adapter)));
}

auto
adapter_matcher::apply_groups(const run_mode &mode) -> void {
  const auto groups =
    get_default_base_groups(std::size(adap_counts), mode.do_groups());
  apply_base_groups(groups, adap_counts);
}
