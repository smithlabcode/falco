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

#include "results_summary.hpp"

#include "falco_utils.hpp"
#include "html.hpp"
#include "report.hpp"
#include "tile_processor.hpp"

#include <algorithm>
#include <array>
#include <assert.h>
#include <iterator>
#include <numeric>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

auto
results_summary::initialize() -> void {
  // assign scalar variables
  const auto gt0 = [](const auto c) { return c > 0; };
  min_read_len = std::ranges::distance(std::cbegin(lengths),
                                       std::ranges::find_if(lengths, gt0));
  total_bases = tabular_dot(lengths);
  const auto gc_acc = [](const auto a, const auto &nuc) {
    return a + nuc[guanine_index] + nuc[cytosine_index];
  };
  total_gc = std::accumulate(std::cbegin(base_counts), std::cend(base_counts),
                             0ul, gc_acc);

  // get summary structures
  centered = tp.get_centered();
  kmer_results = kc.get_kmer_results();

  // apply groups
  groups = get_default_base_groups(max_read_len, do_groups(mode));
  if (do_groups(mode)) {
    apply_base_groups(groups, base_counts);
    apply_base_groups(groups, n_counts);
    // ADS: need one for lengths
    apply_base_groups(groups, qual_by_pos);
  }

  // assign grades
  assign_grades();
}

auto
results_summary::assign_grades() -> void {
  grades.emplace("basic_stats", get_grade_basic_stats());
  grades.emplace("quality_base", get_grade_quality_base(qual_by_pos));
  grades.emplace("quality_sequence", get_grade_quality_sequence(qual_by_read));
  grades.emplace("sequence", get_grade_sequence(base_counts));
  grades.emplace("gc_sequence", get_grade_gc_sequence(gc_content));
  grades.emplace("n_content", get_grade_n_content(n_counts, base_counts));
  grades.emplace("sequence_length", get_grade_sequence_length(lengths));
  grades.emplace("duplication", dr.get_grade_duplication());
  grades.emplace("overrepresented", dr.get_grade_overrepresented());
  grades.emplace("adapter", am.get_grade(n_reads));
  if (do_tiles(mode))
    grades.emplace("tile", get_grade_tile(centered));
  if (do_kmers(mode))
    grades.emplace("kmer", get_grade_kmer(kmer_results));
}

[[nodiscard]] auto
results_summary::get_report() const -> std::string {
  const auto basic_stats = basic_stats_report(
    info, n_reads, min_read_len, max_read_len, total_gc, total_bases, grades);
  auto sections = std::unordered_map<std::string, std::string>{
    {"basic_stats", basic_stats},
    {"quality_base", quality_base_report(qual_by_pos, groups, grades)},
    {"quality_sequence", quality_sequence_report(qual_by_read, grades)},
    {"sequence", sequence_report(base_counts, groups, grades)},
    {"gc_sequence", gc_sequence_report(gc_content, grades)},
    {"n_content", n_content_report(n_counts, base_counts, groups, grades)},
    {"sequence_length", sequence_length_report(lengths, grades)},
    {"duplication", dr.duplication_report(grades)},
    {"overrepresented", dr.overrepresented_report(grades)},
    {"adapter", am.report(n_reads, groups, grades)},
  };
  if (do_tiles(mode)) {
    assert(!centered.empty());
    sections.emplace("tile", tile_report(centered, groups, grades));
  }
  if (do_kmers(mode)) {
    assert(!kmer_results.empty());
    sections.emplace("kmer", kmer_report(kmer_results, grades));
  }
  std::string report;
  for (const auto &name : section_names)
    if (const auto itr = sections.find(name); itr != std::cend(sections))
      report += itr->second;
  return report;
}

auto
results_summary::apply_groups(
  [[maybe_unused]] const std::vector<base_group_t> &groups) -> void {
  return;
}

[[nodiscard]] auto
results_summary::get_html() const -> std::string {
  const auto basic_stats = basic_stats_html(
    info, n_reads, min_read_len, max_read_len, total_gc, total_bases, grades);
  auto sections = std::unordered_map<std::string, std::string>{
    {"basic_stats", basic_stats},
    {"quality_base", quality_base_html(qual_by_pos, groups, grades)},
    {"quality_sequence", quality_sequence_html(qual_by_read, grades)},
    {"sequence", sequence_html(base_counts, groups, grades)},
    {"gc_sequence", gc_sequence_html(gc_content, grades)},
    {"n_content", n_content_html(n_counts, base_counts, groups, grades)},
    {"sequence_length", sequence_length_html(lengths, grades)},
    {"duplication", dr.duplication_html(grades)},
    {"overrepresented", dr.overrepresented_html(grades)},
    {"adapter", am.html(n_reads, groups, grades)},
  };
  if (do_tiles(mode)) {
    assert(!centered.empty());
    sections.emplace("tile", tile_html(centered, groups, grades));
  }
  if (do_kmers(mode)) {
    assert(!kmer_results.empty());
    sections.emplace("kmer", kmer_html(kmer_results, grades));
  }

  assert(std::ranges::all_of(
    sections,  // make sure all grades are configured
    [&](const auto &s) { return grades.is_configured(s); },
    [&](const auto &s) { return s.first; }));

  std::string modules;
  for (const auto &name : section_names)
    if (const auto itr = sections.find(name); itr != std::cend(sections))
      modules += itr->second;
  return falco_get_html(info, grades, modules);
}

[[nodiscard]] auto
results_summary::get_summary() const -> std::string {
  return grades.to_string(info.name);
}
