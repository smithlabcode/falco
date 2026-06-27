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

#include "falco_config.hpp"

#include "falco_grade.hpp"
#include "run_mode.hpp"

#include "boost/boost_unordered.hpp"
#include "nlohmann/json.hpp"

#include <algorithm>
#include <array>
#include <fstream>
#include <iterator>
#include <map>
#include <ranges>
#include <stdexcept>
#include <string>
#include <utility>

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

auto
load_config_and_set_graders(const std::string &filename,
                            run_mode &mode) -> void {
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

  // modes
  for (const auto &label : run_mode::labels)
    if (json_in.contains(label) && !json_in[label].contains("ignore"))
      throw std::runtime_error("invalid config: " + to_string(json_in[label]));

  boost::unordered_flat_map<std::string, bool> modes_in;
  try {
    static constexpr auto ignore = "ignore";
    for (const auto &label : run_mode::labels)
      if (json_in.contains(label) && json_in[label].contains(ignore)) {
        if (json_in[label][ignore] != "0" && json_in[label][ignore] != "1")
          throw std::runtime_error("malformed config entry: " +
                                   to_string(json_in[label]));
        modes_in.emplace(label, json_in[label][ignore] == "0");
      }
  }
  catch (const nlohmann::json::exception &e) {
    throw std::runtime_error(
      std::format("failed parsing config file:\n{}", e.what()));
  }
  mode.assign(modes_in);

  // grades
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
      const auto get_cutoff = [&](const std::string &label0) {
        const auto x = json_in[label][label0].get<std::string>();
        const auto beg = std::data(x);
        const auto end = beg + std::size(x);  // NOLINT(*-pointer-arithmetic)
        double val{};
        const auto res = std::from_chars(beg, end, val);
        if (res.ec != std::errc{})
          throw std::runtime_error("error parsing cutoff: " + x);
        return val;
      };
      const auto fail_val = get_cutoff("error");
      const auto warn_val = get_cutoff("warn");
      graders.emplace(label, grader(label, warn_val, fail_val));
    }
  }
  catch (const nlohmann::json::exception &e) {
    throw std::runtime_error(
      std::format("failed parsing config file:\n{}", e.what()));
  }
  grader_set::instance(graders);
}
