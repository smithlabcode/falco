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
#include <ranges>
#include <stdexcept>
#include <string>

[[nodiscard]] auto
analysis_grades::summary(const std::string &infile_path) const -> std::string {
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
  for (const auto [name, label] : std::views::zip(names, labels))
    if (is_configured(name))
      r += compose(grade(name), label);
  return r;
}

[[nodiscard]] auto
analysis_grades::grade(const std::string &name) const -> std::string {
  const auto itr = g.find(name);
  return itr == std::cend(g) ? std::string{} : itr->second;
}

[[nodiscard]] auto
analysis_grades::label(const std::string &name) const -> std::string {
  const auto itr = std::ranges::find(names, name);
  if (itr == std::cend(names))
    throw std::runtime_error(std::format("bad grade category: {}", name));
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
  return labels[std::distance(std::cbegin(names), itr)];
}
