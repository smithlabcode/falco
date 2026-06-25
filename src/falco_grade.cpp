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

#include <cctype>
#include <cstdint>
#include <filesystem>
#include <format>
#include <fstream>
#include <ranges>
#include <stdexcept>
#include <string>
#include <tuple>

// clang-format off
static const auto default_graders = std::array{
  grader{.name="quality_base_median", .warn=25.0, .fail=20.0},
  grader{.name="quality_base_lower",  .warn=10.0, .fail=5.00},
  grader{.name="quality_sequence",    .warn=27.0, .fail=20.0},
  grader{.name="n_content",           .warn=0.05, .fail=0.20},
  grader{.name="sequence",            .warn=10.0, .fail=20.0},
  grader{.name="sequence_length",     .warn=1.00, .fail=1.00},
  grader{.name="gc_sequence",         .warn=15.0, .fail=30.0},
  grader{.name="duplication",         .warn=0.70, .fail=0.50},
  grader{.name="overrepresented",     .warn=0.10, .fail=1.00},
  grader{.name="kmer",                .warn=2.00, .fail=5.00},
  grader{.name="tile",                .warn=5.00, .fail=10.0},
  grader{.name="adapter",             .warn=0.05, .fail=0.10},
};
// clang-format on

[[nodiscard]] auto
file_grades::summary(const std::string &infile_path) const -> std::string {
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

[[nodiscard]] auto
strip(const std::string &s) -> std::string {
  if (s.empty())
    return {};
  const auto is_print = [](const auto c) { return std::isprint(c); };
  const auto start = std::ranges::find_if(s, is_print);
  if (start == std::cend(s))
    return std::string{};
  const auto stop = std::ranges::find_last_if(s, is_print);
  return std::string{start, std::cbegin(stop) + 1};
}

[[nodiscard]] auto
split(const std::string &s) {
  std::istringstream iss(s);
  return std::vector<std::string>{std::istream_iterator<std::string>{iss}, {}};
}

[[nodiscard]] static auto
load_grades(const std::string &filename)
  -> boost::unordered_flat_map<std::string, grader> {
  // ADS: (todo) handle carriage returns and other control chars
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("failed to open config file: " + filename);
  nlohmann::json json_in;

  std::string line;
  while (std::getline(in, line)) {
    line = strip(line);
    if (line.empty() || line[0] == '#')
      continue;
    const auto parts = split(line);
    if (std::size(parts) != 3)
      throw std::runtime_error("malformed config line: " + line);
    json_in[parts[0]][parts[1]] = parts[2];
  }

  for (const auto &label : grade_labels)
    if (json_in.contains(label) &&
        !(json_in[label].contains("error") && json_in[label].contains("warn")))
      throw std::runtime_error("missing config value for: " +
                               std::string(label));

  boost::unordered_flat_map<std::string, grader> graders;
  try {
    for (const auto &label : grade_labels) {
      if (!json_in.contains(label))
        continue;
      const auto fail_elem = json_in[label]["error"].get<std::string>();
      const auto warn_elem = json_in[label]["warn"].get<std::string>();
      // ADS: fix this below
      const auto fail_val = std::atof(std::data(fail_elem));
      const auto warn_val = std::atof(std::data(warn_elem));
      graders.emplace(label, grader(label, warn_val, fail_val));
    }
  }
  catch (const nlohmann::json::exception &e) {
    throw std::runtime_error(
      std::format("failed parsing config file:\n{}", e.what()));
  }
  return graders;
}

grader_set::grader_set(const std::string &filename) {
  if (filename.empty()) {
    using map_t = boost::unordered_flat_map<std::string, grader>;
    const auto make_elem = [](const auto &x) { return std::pair{x.name, x}; };
    graders = default_graders | std::views::transform(make_elem) |
              std::ranges::to<map_t>();
  }
  else
    graders = load_grades(filename);
}
