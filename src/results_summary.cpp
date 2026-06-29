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
results_summary::apply_groups() -> void {
  groups = get_default_base_groups(max_read_len, mode.do_groups());
  if (mode.do_groups()) {
    apply_base_groups(groups, base_counts);
    apply_base_groups(groups, n_counts);
    // ADS: need one for lengths
    apply_base_groups(groups, qual_by_pos);
    am.apply_groups(mode);
    tp.apply_groups(mode);
  }
}

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
  median_read_len = median_tabular(lengths);

  // apply groups before making summary stats like in FastQC
  apply_groups();

  // get summary structures
#ifdef ORIGINAL_DUPS
  n_reads_for_dups = n_reads;
  dup_summary = dr.get_dups_summary(n_reads_for_dups);
  dup_summary.n_reads = n_reads_for_dups;
#else   //  ORIGINAL_DUPS
  n_reads_for_dups = dr.get_n_reads();
  dup_summary = dr.get_dups_summary();
#endif  //  ORIGINAL_DUPS
  overrep = dr.get_overrepresented(n_reads_for_dups);
  centered = tp.get_centered();
  kmer_results = kc.get_kmer_results();

  // assign grades
  assign_grades();
}

auto
results_summary::assign_grades() -> void {
  grades.emplace("basic_stats", get_grade_basic_stats());

  if (mode.do_adap())
    grades.emplace("adapter", am.get_grade(n_reads));

  if (mode.do_dups())
    grades.emplace("duplication", get_grade_duplication(dup_summary));

  if (mode.do_gc_content())
    grades.emplace("gc_sequence", get_grade_gc_sequence(gc_content));

  if (mode.do_kmers())
    grades.emplace("kmer", get_grade_kmer(kmer_results));

  if (mode.do_length())
    grades.emplace("sequence_length", get_grade_sequence_length(lengths));

  if (mode.do_n_content())
    grades.emplace("n_content", get_grade_n_content(n_counts, base_counts));

  if (mode.do_overrep())
    grades.emplace("overrepresented",
                   get_grade_overrepresented(n_reads_for_dups, dr));

  if (mode.do_qual_base())
    grades.emplace("quality_base", get_grade_quality_base(qual_by_pos));

  if (mode.do_qual_seq())
    grades.emplace("quality_sequence",
                   get_grade_quality_sequence(qual_by_read));

  if (mode.do_sequence())
    grades.emplace("sequence", get_grade_sequence(base_counts));

  if (mode.do_tiles())
    grades.emplace("tile", get_grade_tile(centered));
}

[[nodiscard]] auto
results_summary::get_report() const -> std::string {
  const auto basic_stats =
    basic_stats_report(info, n_reads, min_read_len, max_read_len,
                       median_read_len, total_gc, total_bases, grades);
  auto sections = std::unordered_map<std::string, std::string>{
    {"basic_stats", basic_stats},
  };

  if (mode.do_adap())
    sections.emplace("adapter", am.report(n_reads, groups, grades));

  if (mode.do_dups())
    sections.emplace("duplication", duplication_report(dup_summary, grades));

  if (mode.do_gc_content())
    sections.emplace("gc_sequence", gc_sequence_report(gc_content, grades));

  if (mode.do_kmers()) {
    assert(!kmer_results.empty());
    sections.emplace("kmer", kmer_report(kmer_results, grades));
  }

  if (mode.do_length())
    sections.emplace("sequence_length",
                     sequence_length_report(lengths, grades));

  if (mode.do_n_content())
    sections.emplace("n_content",
                     n_content_report(n_counts, base_counts, groups, grades));

  if (mode.do_overrep())
    sections.emplace("overrepresented",
                     overrepresented_report(overrep, grades));

  if (mode.do_qual_base())
    sections.emplace("quality_base",
                     quality_base_report(qual_by_pos, groups, grades));

  if (mode.do_qual_seq())
    sections.emplace("quality_sequence",
                     quality_sequence_report(qual_by_read, grades));

  if (mode.do_sequence())
    sections.emplace("sequence", sequence_report(base_counts, groups, grades));

  if (mode.do_tiles()) {
    assert(!centered.empty());
    sections.emplace("tile", tile_report(centered, groups, grades));
  }

  std::string report;
  for (const auto &name : section_names)
    if (const auto itr = sections.find(name); itr != std::cend(sections))
      report += itr->second;
  return report;
}

[[nodiscard]] auto
results_summary::get_html() const -> std::string {
  const auto basic_stats = basic_stats_html(
    info, n_reads, min_read_len, max_read_len, total_gc, total_bases, grades);
  auto sections = std::unordered_map<std::string, std::string>{
    {"basic_stats", basic_stats},
  };

  if (mode.do_adap())
    sections.emplace("adapter", am.html(n_reads, groups, grades));

  if (mode.do_dups())
    sections.emplace("duplication", duplication_html(dup_summary, grades));

  if (mode.do_gc_content())
    sections.emplace("gc_sequence", gc_sequence_html(gc_content, grades));

  if (mode.do_kmers()) {
    assert(!kmer_results.empty());
    sections.emplace("kmer", kmer_html(kmer_results, grades));
  }

  if (mode.do_length())
    sections.emplace("sequence_length", sequence_length_html(lengths, grades));

  if (mode.do_n_content())
    sections.emplace("n_content",
                     n_content_html(n_counts, base_counts, groups, grades));

  if (mode.do_overrep())
    sections.emplace("overrepresented", overrepresented_html(overrep, grades));

  if (mode.do_qual_base())
    sections.emplace("quality_base",
                     quality_base_html(qual_by_pos, groups, grades));

  if (mode.do_qual_seq())
    sections.emplace("quality_sequence",
                     quality_sequence_html(qual_by_read, grades));

  if (mode.do_sequence())
    sections.emplace("sequence", sequence_html(base_counts, groups, grades));

  if (mode.do_tiles()) {
    assert(!centered.empty());
    sections.emplace("tile", tile_html(centered, groups, grades));
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
