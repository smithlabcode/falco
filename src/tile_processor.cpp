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

#include "tile_processor.hpp"
#include "falco_utils.hpp"

#include <algorithm>
#include <cstdint>
#include <fstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

std::uint32_t tile_processor::preceding_colons = 0;

auto
tile_processor::set_preceding_colons(const std::string &fastq_buffername)
  -> std::uint32_t {
  std::ifstream in(fastq_buffername);
  if (!in)
    throw std::runtime_error("failed to open file: " + fastq_buffername);
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("failed to read line from: " + fastq_buffername);
  const auto colons_found = std::ranges::count(line, ':');
  preceding_colons = colons_found >= 6 ? 4 : (colons_found >= 4 ? 2 : 0);
  return preceding_colons;
}

[[nodiscard]] auto
tile_processor::string() const -> std::string {
  static constexpr auto as_frac = [](const auto a, const auto b) {
    return static_cast<double>(a) / static_cast<double>(b);
  };
  auto r = std::format(">>Per tile sequence quality\t{}\n", "pass");
  r += header;
  for (const auto &[i, q] : quals)
    for (auto j = 0; j < std::size(q); ++j)
      r += std::format("{}\t{}\t{:.6g}\n", i, j + 1,
                       as_frac(q[j].first, q[j].second));
  r += footer;
  return r + end_module_tag;
}

auto
tile_processor::operator+=(const tile_processor &rhs)
  -> const tile_processor & {
  const auto pair_plus = [](const auto &a, const auto &b) {
    return std::pair{a.first + b.first, a.second + b.second};
  };
  for (const auto &[i, q] : rhs.quals) {
    if (quals.contains(i)) {
      quals[i].resize(std::max(std::size(quals[i]), std::size(q)));
      std::ranges::transform(quals[i], q, std::begin(quals[i]), pair_plus);
    }
    else
      quals.emplace(i, q);
  }
  return *this;
}
