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

#include "falco_grade.hpp"

#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <limits>
#include <numeric>
#include <ranges>
#include <sstream>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <vector>

// clang-format off
static const auto default_graders = std::array{  // NOLINT(cert-err58-cpp)
  grader{.name = "quality_base_median", .warn = 25.0, .fail = 20.0},
  grader{.name = "quality_base_lower",  .warn = 10.0, .fail = 5.00},
  grader{.name = "quality_sequence",    .warn = 27.0, .fail = 20.0},
  grader{.name = "n_content",           .warn = 0.05, .fail = 0.20},
  grader{.name = "sequence",            .warn = 10.0, .fail = 20.0},
  grader{.name = "sequence_length",     .warn = 1.00, .fail = 1.00},
  grader{.name = "gc_sequence",         .warn = 15.0, .fail = 30.0},
  grader{.name = "duplication",         .warn = 0.70, .fail = 0.50},
  grader{.name = "overrepresented",     .warn = 0.10, .fail = 1.00},
  grader{.name = "kmer",                .warn = 2.00, .fail = 5.00},
  grader{.name = "tile",                .warn = 5.00, .fail = 10.0},
  grader{.name = "adapter",             .warn = 0.05, .fail = 0.10},
};
// clang-format on

[[nodiscard]] auto
file_grades::to_string(const std::string &infile_path) const -> std::string {
  const auto infile = std::filesystem::path(infile_path).filename().string();
  const auto up = [](const auto &s) {
    const auto u = [](const std::uint8_t c) { return std::toupper(c); };
    return std::views::transform(s, u) | std::ranges::to<std::string>();
  };
  const auto compose = [&](const auto &g, const auto &lbl) {
    static constexpr auto f = "{}\t{}\t{}\n";
    return g.empty() ? std::string{} : std::format(f, up(g), lbl, infile);
  };
  std::string r;
  for (const auto [name, title] :
       std::views::zip(section_names, section_titles))
    if (is_configured(name))
      r += compose(grade(name), title);
  return r;
}

[[nodiscard]] auto
file_grades::grade(const std::string &label) const -> std::string {
  const auto itr = g.find(label);
  return itr == std::cend(g) ? std::string{} : itr->second;
}

[[nodiscard]] auto
file_grades::get_title(const std::string &name) const -> std::string {
  const auto itr = std::ranges::find(section_names, name);
  if (itr == std::cend(section_names))
    throw std::runtime_error(std::format("bad section name: {}", name));
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  return section_titles[std::distance(std::cbegin(section_names), itr)];
}

[[nodiscard]] auto
grader::identify_grade(const double value) const -> std::string {
  if (warn < fail) {
    if (value < warn)
      return "pass";
    if (value < fail)
      return "warn";
    return "fail";
  }
  else {
    if (value < fail)
      return "fail";
    if (value < warn)
      return "warn";
    return "pass";
  }
  throw std::runtime_error(
    std::format("value {} outside grade range: {}", value, *this));
}

[[nodiscard]] auto
grader::to_string() const -> std::string {
  static constexpr auto n_indent = 4;
  nlohmann::json data = *this;
  return data.dump(n_indent);
}

[[nodiscard]] auto
grader_set::get_grader(const std::string &label) -> const grader & {
  const auto &g = instance().graders;
  const auto itr = g.find(label);
  if (itr == std::cend(g))
    throw std::runtime_error("missing grader for: " + label);
  return itr->second;
}

[[nodiscard]] auto
grader_set::get_grade(const std::string &label,
                      const double value) -> std::string {
  return get_grader(label).identify_grade(value);
}

grader_set::grader_set(const map_t<std::string, grader> &g) {
  if (g.empty()) {
    const auto make_elem = [](const auto &x) { return std::pair{x.name, x}; };
    graders = default_graders | std::views::transform(make_elem) |
              std::ranges::to<map_t<std::string, grader>>();
  }
  else
    graders = g;
}

[[nodiscard]] auto
get_grade_sequence_length(const std::vector<std::uint64_t> &lengths)
  -> std::string {
  // ADS: this module just has binary toggles which need to be
  // incorporated
  [[maybe_unused]] static constexpr auto label = "sequence_length";
  const bool has_empty_reads = std::size(lengths) > 0 && lengths[0] > 0;
  if (has_empty_reads)
    return "fail";
  const auto n_lengths =
    std::ranges::count_if(lengths, [](const auto x) { return x > 0; });
  if (n_lengths > 1)
    return "warn";
  return "pass";
}

[[nodiscard]] auto
get_grade_gc_sequence(const falco::gc_content_array &gc_content)
  -> std::string {
  static constexpr auto label = "gc_sequence";
  return grader_set::get_grade(label, sum_deviation_from_normal(gc_content));
}

[[nodiscard]] auto
get_grade_sequence(const std::vector<falco::nuc_array> &nucs) -> std::string {
  static constexpr auto label = "sequence";
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
  return grader_set::get_grade(label, max_diff);
}

[[nodiscard]] auto
get_grade_n_content(const std::vector<std::uint64_t> &n_counts,
                    const std::vector<falco::nuc_array> &nucs) -> std::string {
  static constexpr auto label = "n_content";
  // ADS: at this point 'nucs' should not include counts of 'N'
  const auto max_idx =
    std::distance(std::cbegin(n_counts), std::ranges::max_element(n_counts));
  const auto total_non_n =
    std::reduce(std::cbegin(nucs[max_idx]), std::cend(nucs[max_idx]));
  const auto max_n_frac =
    as_frac(n_counts[max_idx], n_counts[max_idx] + total_non_n);
  return grader_set::get_grade(label, max_n_frac);
}

[[nodiscard]] auto
get_grade_quality_sequence(const falco::qual_array &qual_by_read)
  -> std::string {
  static constexpr auto label = "quality_sequence";
  const auto q_beg = std::cbegin(qual_by_read);
  const auto max_itr = std::ranges::max_element(qual_by_read);
  const auto qual_val_mode = std::distance(q_beg, max_itr);
  return grader_set::get_grade(label, static_cast<double>(qual_val_mode));
}

[[nodiscard]] auto
get_grade_quality_base(const std::vector<falco::qual_array> &qual)
  -> std::string {
  static constexpr auto label_median = "quality_base_median";
  static constexpr auto label_lquart = "quality_base_lower";
  auto min_qual_median = std::numeric_limits<std::uint32_t>::max();
  auto min_qual_lquart = std::numeric_limits<std::uint32_t>::max();
  for (const auto &q : qual) {
    const auto quantiles = five_quants(q);
    min_qual_median = std::min(min_qual_median, median_val(quantiles));
    min_qual_lquart = std::min(min_qual_lquart, lquart_val(quantiles));
  }
  const auto lq_grade = grader_set::get_grade(label_lquart, min_qual_lquart);
  const auto med_grade = grader_set::get_grade(label_median, min_qual_median);
  using std::string_literals::operator""s;
  return (lq_grade == "fail"s || med_grade == "fail"s)   ? "fail"s
         : (lq_grade == "warn"s || med_grade == "warn"s) ? "warn"s
                                                         : "pass"s;
}

[[nodiscard]] auto
get_grade_basic_stats() -> std::string {
  static constexpr auto default_grade = "pass";  // always a pass
  return default_grade;
}
