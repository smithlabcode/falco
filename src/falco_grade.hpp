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

#include "falco_utils.hpp"
#include "quality_score.hpp"

#include "boost/boost_unordered.hpp"
#include "nlohmann/json.hpp"

#include <array>
#include <cstdint>
#include <format>
#include <iterator>
#include <string>
#include <vector>  // IWYU pragma: keep

// clang-format off

static constexpr auto section_names = std::array{
  "basic_stats",
  "quality_base",
  "tile",
  "quality_sequence",
  "sequence",
  "gc_sequence",
  "n_content",
  "sequence_length",
  "duplication",
  "overrepresented",
  "adapter",
  "kmer",
};

static constexpr auto section_titles = std::array{
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

static constexpr auto grade_labels = std::array{
  "basic_stats",
  "quality_base_median",
  "quality_base_lower",
  "quality_sequence",
  "tile",
  "sequence",
  "gc_sequence",
  "n_content",
  "sequence_length",
  "duplication",
  "overrepresented",
  "adapter",
  "kmer",
};

static_assert(std::size(section_names) == std::size(section_titles));
// clang-format on

struct file_grades {
  boost::unordered_flat_map<std::string, std::string> g;

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
  get_title(const std::string &name) const -> std::string;

  [[nodiscard]] auto
  to_string(const std::string &infile_path) const -> std::string;
};

struct grader {
  std::string name;
  double warn{};
  double fail{};

  [[nodiscard]] auto
  identify_grade(const double value) const -> std::string;

  [[nodiscard]] auto
  to_string() const -> std::string;
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(grader, name, warn, fail);
};

template <> struct std::formatter<grader> : std::formatter<std::string> {
  auto
  format(const grader &g, auto &ctx) const {
    return std::formatter<std::string>::format(g.to_string(), ctx);
  }
};

class grader_set {
  template <typename T, typename U>
  using map_t = boost::unordered_flat_map<T, U>;

public:
  static auto
  instance(const map_t<std::string, grader> &g = {}) -> const grader_set & {
    static const grader_set s(g);
    return s;
  }

  [[nodiscard]] static auto
  get_grader(const std::string &label) -> const grader &;

  [[nodiscard]] static auto
  get_grade(const std::string &label, const double value) -> std::string;

private:
  explicit grader_set(const map_t<std::string, grader> &g);
  grader_set() = default;
  ~grader_set() = default;

  map_t<std::string, grader> graders;
};  // grader_set

[[nodiscard]] auto
get_grade_sequence_length(const std::vector<std::uint64_t> &lengths)
  -> std::string;

[[nodiscard]] auto
get_grade_gc_sequence(const falco::gc_content_array &gc_content) -> std::string;

[[nodiscard]] auto
get_grade_sequence(const std::vector<falco::nuc_array> &nucs) -> std::string;

[[nodiscard]] auto
get_grade_n_content(const std::vector<std::uint64_t> &n_counts,
                    const std::vector<falco::nuc_array> &nucs) -> std::string;

[[nodiscard]] auto
get_grade_quality_sequence(const falco::qual_array &qual_by_read)
  -> std::string;

[[nodiscard]] auto
get_grade_quality_base(const std::vector<falco::qual_array> &qual)
  -> std::string;

[[nodiscard]] auto
get_grade_basic_stats() -> std::string;

#endif  // SRC_FALCO_GRADE_HPP_
