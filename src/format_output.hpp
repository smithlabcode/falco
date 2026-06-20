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

#ifndef SRC_FORMAT_OUTPUT_HPP_
#define SRC_FORMAT_OUTPUT_HPP_

#include "falco_utils.hpp"
#include "quality_score.hpp"

#include <cstdint>
#include <format>
#include <limits>  // IWYU pragma: keep
#include <string>
#include <vector>  // IWYU pragma: keep

// clang-format off
static constexpr auto report_section_order = std::array{
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
  "kmers",
};
// clang-format on

[[nodiscard]] auto
get_grade_read_lengths(const std::vector<std::uint64_t> &lengths)
  -> std::string;

[[nodiscard]] auto
format_read_lengths(const std::vector<std::uint64_t> &lengths,
                    const std::string &grade) -> std::string;

[[nodiscard]] auto
format_read_lengths_html(const std::vector<std::uint64_t> &lengths,
                         const std::string &grade) -> std::string;

[[nodiscard]] auto
get_grade_gc_content(const falco::gc_content_array &gc_content) -> std::string;

[[nodiscard]] auto
format_gc_content(const falco::gc_content_array &gc_content,
                  const std::string &grade) -> std::string;

[[nodiscard]] auto
format_gc_content_html(const falco::gc_content_array &gc_content,
                       const std::string &grade) -> std::string;

[[nodiscard]] auto
get_grade_base_comp(const std::vector<falco::nuc_array> &nucs) -> std::string;

[[nodiscard]] auto
format_base_comp(const std::vector<falco::nuc_array> &nucs,
                 const std::vector<base_group_t> &groups,
                 const std::string &grade) -> std::string;

[[nodiscard]] auto
format_base_comp_html(const std::vector<falco::nuc_array> &nucs,
                      const std::vector<base_group_t> &groups,
                      const std::string &grade) -> std::string;

[[nodiscard]] auto
get_grade_n_content(const std::vector<std::uint64_t> &n_counts,
                    const std::vector<falco::nuc_array> &nucs) -> std::string;

[[nodiscard]] auto
format_n_content(const std::vector<std::uint64_t> &n_counts,
                 const std::vector<falco::nuc_array> &nucs,
                 const std::vector<base_group_t> &groups,
                 const std::string &grade) -> std::string;

[[nodiscard]] auto
format_n_content_html(const std::vector<std::uint64_t> &n_counts,
                      const std::vector<falco::nuc_array> &nucs,
                      const std::vector<base_group_t> &groups,
                      const std::string &grade) -> std::string;

[[nodiscard]] auto
get_grade_qual_by_read(const falco::qual_array &qual_by_read) -> std::string;

[[nodiscard]] auto
format_qual_by_read(const falco::qual_array &qual_by_read,
                    const std::string &grade) -> std::string;

[[nodiscard]] auto
format_qual_by_read_html(const falco::qual_array &qual_by_read,
                         const std::string &grade) -> std::string;

[[nodiscard]] auto
get_grade_qual_by_pos(const std::vector<falco::qual_array> &qual)
  -> std::string;

[[nodiscard]] auto
format_qual_by_pos(const std::vector<falco::qual_array> &qual,
                   const std::vector<base_group_t> &groups,
                   const std::string &grade) -> std::string;

[[nodiscard]] auto
format_qual_by_pos_html(const std::vector<falco::qual_array> &qual,
                        const std::vector<base_group_t> &groups,
                        const std::string &grade) -> std::string;

[[nodiscard]] auto
get_grade_basic_stats() -> std::string;

[[nodiscard]] auto
format_basic_stats(const file_info &info, const std::uint64_t n_reads,
                   const std::uint64_t min_read_len,
                   const std::uint64_t max_read_len,
                   const std::uint64_t total_gc, const std::uint64_t total_nucs,
                   const std::string &grade) -> std::string;

[[nodiscard]] auto
format_basic_stats_html(const file_info &info, const std::uint64_t n_reads,
                        const std::uint64_t min_read_len,
                        const std::uint64_t max_read_len,
                        const std::uint64_t total_gc,
                        const std::uint64_t total_nucs) -> std::string;

#endif  // SRC_FORMAT_OUTPUT_HPP_
