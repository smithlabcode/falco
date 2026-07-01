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

#ifndef SRC_HTML_HPP_
#define SRC_HTML_HPP_

#include "falco_utils.hpp"
#include "quality_score.hpp"
#include "tile_processor.hpp"

#include <cstdint>
#include <format>
#include <string>

struct file_grades;
struct kmer_result;

// grade / section_label / section_title / grade / text
static constexpr auto html_module_fmt =
  R"(<div class="module">
<h2 class="{}" id="{}">{}: {}</h2>
{}
</div>)";

[[nodiscard]] auto
falco_get_html(const file_info &info, const file_grades &grades,
               const std::string &analysis_modules) -> std::string;

[[nodiscard]] auto
sequence_length_html(const std::vector<std::uint64_t> &lengths,
                     const file_grades &grades) -> std::string;

[[nodiscard]] auto
gc_sequence_html(const falco::gc_content_array &gc_content,
                 const file_grades &grades) -> std::string;

[[nodiscard]] auto
sequence_html(const std::vector<falco::nuc_array> &nucs,
              const std::vector<base_group_t> &groups,
              const file_grades &grades) -> std::string;

[[nodiscard]] auto
n_content_html(const std::vector<std::uint64_t> &n_counts,
               const std::vector<falco::nuc_array> &nucs,
               const std::vector<base_group_t> &groups,
               const file_grades &grades) -> std::string;

[[nodiscard]] auto
quality_sequence_html(const falco::qual_array &qual_by_read,
                      const file_grades &grades) -> std::string;

[[nodiscard]] auto
quality_base_html(const std::vector<falco::qual_array> &qual,
                  const std::vector<base_group_t> &groups,
                  const file_grades &grades) -> std::string;

[[nodiscard]] auto
basic_stats_html(const file_info &info, const std::uint64_t n_reads,
                 const std::uint64_t min_read_len,
                 const std::uint64_t max_read_len,
                 const std::uint64_t median_read_len,
                 const std::uint64_t total_gc, const std::uint64_t total_nucs,
                 const file_grades &grades) -> std::string;

[[nodiscard]] auto
tile_html(const tile_processor::tiles_centered_t &centered,
          const std::vector<base_group_t> &groups,
          const file_grades &grades) -> std::string;

[[nodiscard]] auto
kmer_html(const std::vector<kmer_result> &results,
          const file_grades &grades) -> std::string;

#endif  // SRC_HTML_HPP_
