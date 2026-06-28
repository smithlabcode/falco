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

#include "adapter_set.hpp"
#include "run_mode.hpp"

#include <algorithm>
#include <cctype>
#include <format>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

[[nodiscard]] static auto
load_adapters(const std::string &filename) {
  // ADS: (todo) handle carriage returns and other control chars
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("failed to open adapters file: " + filename);

  std::vector<std::string> adapter_names;
  std::vector<std::string> adapters;

  std::string line_data;
  while (std::getline(in, line_data)) {
    std::string_view line = line_data;
    const auto to_keep_prefix = line.find_first_not_of(" \t");
    if (to_keep_prefix == std::string_view::npos)
      continue;
    line.remove_prefix(std::min(to_keep_prefix, std::size(line)));
    if (line[0] == '#')
      continue;
    const auto to_keep_suffix = line.find_last_not_of(" \t");
    if (to_keep_suffix == std::string_view::npos)
      continue;
    line.remove_suffix(std::size(line) - to_keep_suffix - 1);
    std::string cleaned_line;
    for (auto itr = std::cbegin(line); itr != std::cend(line); ++itr)
      if (!std::isblank(*itr) ||
          (std::next(itr) != std::cend(line) && *itr != *std::next(itr)))
        cleaned_line += *itr;
    const auto tab_pos = cleaned_line.find('\t');
    if (tab_pos == std::string::npos ||
        tab_pos != cleaned_line.find_last_of('\t'))
      throw std::runtime_error("malformed line: " + line_data);
    const auto is_print = [](const auto c) { return std::isprint(c); };
    const auto name = cleaned_line.substr(0, tab_pos);
    const auto seq = cleaned_line.substr(tab_pos + 1);
    if (!std::ranges::all_of(name, is_print) ||
        !std::ranges::all_of(seq, is_print))
      throw std::runtime_error("malformed line: " + line_data);

    adapter_names.push_back(name);
    adapters.push_back(seq);
  }
  return std::tuple{adapter_names, adapters};
}

adapter_set::adapter_set(const run_mode &m, const std::string &filename) {
  if (m.do_adap()) {  // ADS: if not do_adap, then this is empty
    if (filename.empty()) {
      std::ranges::copy(default_adapter_names,
                        std::back_inserter(adapter_names));
      std::ranges::copy(default_adapters, std::back_inserter(adapters));
    }
    else
      std::tie(adapter_names, adapters) = load_adapters(filename);
  }
}

[[nodiscard]] auto
adapter_set::validate() const -> std::tuple<bool, std::string> {
  const auto adap_len = std::size(instance().adapters.front());
  if (!std::ranges::all_of(instance().adapters, [&](const auto &a) {
        return std::size(a) == adap_len;
      }))
    return {false, "adapter patterns must have equal length"};
  if (adap_len < min_adapter_size)
    return {false, std::format("adapter patterns must have length at least {}",
                               min_adapter_size)};
  if (adap_len > max_adapter_size)
    return {false, std::format("adapter patterns must have length at most {}",
                               max_adapter_size)};
  return {true, std::string{}};
}
