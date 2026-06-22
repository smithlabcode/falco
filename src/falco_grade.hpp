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

#ifndef SRC_FALCO_GRADE_HPP_
#define SRC_FALCO_GRADE_HPP_

#include <algorithm>
#include <array>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

struct analysis_grades {
  // clang-format off
  static constexpr auto names = std::array{
    "basic_stats",
    "qual_by_pos",
    "tiles",
    "qual_by_read",
    "base_comp",
    "gc_content",
    "n_content",
    "read_lengths",
    "duplication",
    "overrep",
    "adapters",
    "kmer",
  };
  static constexpr auto labels = std::array{
    "Basic Statistics",
    "Per base sequence quality",
    "Per tile sequence quality",
    "Per sequence quality scores",
    "Per base sequence content",
    "Per sequence GC content",
    "Per base N content",
    "Sequence Length Distribution",
    "Sequence Duplication Levels",
    "Overrepresented sequences",
    "Adapter Content",
    "Kmer Content",
  };
  static_assert(std::size(names) == std::size(labels));
  // clang-format on
  std::unordered_map<std::string, std::string> g;

  [[nodiscard]] static auto
  is_valid(const std::string &name) -> bool {
    return std::ranges::find(labels, name) != std::cend(labels);
  }

  [[nodiscard]] auto
  is_configured(const std::string &name) const -> bool {
    return g.contains(name);
  }

  auto
  emplace(const std::string &name, const std::string &grade) {
    g.emplace(name, grade);
  }

  [[nodiscard]] auto
  grade(const std::string &name) const -> std::string;

  [[nodiscard]] auto
  get_label(const std::string &name) const -> std::string;

  [[nodiscard]] auto
  summary(const std::string &infile_path) const -> std::string;
};

[[nodiscard]] inline auto
identify_grade(const auto &cutoffs, const auto c) {
  const auto a = std::pair{c, std::string{}};
  const auto b = std::ranges::lower_bound(cutoffs, a);
  if (b == std::cend(cutoffs))
    throw std::runtime_error("error in identifying grade");
  return b->second;
}

#endif  // SRC_FALCO_GRADE_HPP_
