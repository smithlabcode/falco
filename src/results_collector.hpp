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

#ifndef SRC_RESULTS_COLLECTOR_HPP_
#define SRC_RESULTS_COLLECTOR_HPP_

#include "adapter_matcher.hpp"
#include "bam_file.hpp"
#include "duplication_results.hpp"
#include "falco_file_format.hpp"
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "format_output.hpp"
#include "html.hpp"
#include "kmer_counter.hpp"
#include "tile_processor.hpp"

#include <algorithm>
#include <cstdint>
#include <ranges>
#include <string>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <vector>

[[nodiscard]] static inline auto
array_contains(const auto &a, const auto &b) {
  return std::ranges::find(a, b) != std::cend(a);
};

static constexpr auto assumed_page_size = 4096;
struct alignas(assumed_page_size) results_collector {
  // struct alignas(std::hardware_destructive_interference_size)
  // results_collector {
  std::uint64_t n_reads{};
  std::uint64_t max_read_len{};
  std::vector<falco::nuc_array> base_counts;
  falco::gc_content_array gc_content{};
  std::vector<std::uint64_t> n_counts;
  std::vector<std::uint64_t> lengths;
  std::vector<falco::qual_array> qual_by_pos;
  falco::qual_array qual_by_read{};
  duplication_results dr;
  adapter_matcher am;
  std::string seq;

  results_collector() : lengths(1, 0) {}  // in case all reads have length 0

  auto
  resize(const std::uint32_t updated_length) {
    base_counts.resize(updated_length);
    n_counts.resize(updated_length);
    lengths.resize(updated_length + 1);  // need one extra here
    qual_by_pos.resize(updated_length);
    am.resize(updated_length);
  }

  // clang-format off
  template <typename self_t>
  auto
  init(this self_t &self, const auto &info) { self.init_impl(info); }
  auto init_impl(const auto &info) { dr.initialize(info.n_reads_est); };
  // clang-format on

  template <typename self_t>
  auto
  get_grades(this self_t &self) -> file_grades {
    file_grades grades;
    self.get_grades_impl(grades);
    return grades;
  }

  auto
  get_grades_impl(file_grades &grades) {
    grades.emplace("basic_stats", get_grade_basic_stats());
    grades.emplace("quality_base", get_grade_quality_base(qual_by_pos));
    grades.emplace("quality_sequence",
                   get_grade_quality_sequence(qual_by_read));
    grades.emplace("sequence", get_grade_sequence(base_counts));
    grades.emplace("gc_sequence", get_grade_gc_sequence(gc_content));
    grades.emplace("n_content", get_grade_n_content(n_counts, base_counts));
    grades.emplace("sequence_length", get_grade_sequence_length(lengths));
    grades.emplace("duplication", dr.get_grade_duplication());
    grades.emplace("overrepresented", dr.get_grade_overrepresented());
    grades.emplace("adapter", am.get_grade(n_reads));
  }

  template <typename self_t>
  auto
  finalize(this self_t &self, const run_mode &mode, const file_info &info) {
    self.finalize_impl(mode, info);
  }

  auto
  finalize_impl(const run_mode &mode, const file_info &info) {
    adjust_base_counts_for_ns();
    if (!is_mapped_reads(info.format))
      adjust_fastq_qual_encoding(qual_by_pos, qual_by_read, info.encoding);
    if (do_groups(mode)) {
      const auto groups =
        get_default_base_groups(max_read_len, do_groups(mode));
      apply_base_groups(groups, base_counts);
      apply_base_groups(groups, n_counts);
      // ADS: need one for lengths
      apply_base_groups(groups, qual_by_pos);
      am.finalize(mode);
    }
  }

  auto
  adjust_base_counts_for_ns()
    -> void {  // Ns were counted among the C in base_counts
    for (auto i = 0u; i < std::size(base_counts); ++i)
      base_counts[i][3] -= n_counts[i];
  }

  template <typename self_t>
  auto
  process_one_read(this self_t &self, const auto &rec) -> void {
    self.process_one_read_impl(rec);
  }

  template <typename rec_t, typename self_t>
  auto
  process_reads(this self_t &self, auto cursor, const auto lim) {
    rec_t rec{};
    while (cursor < lim && (rec = get_next(cursor, lim))) {
      self.process_one_read(rec);
      ++self.n_reads;
    }
  }

  // clang-format off
  [[nodiscard]] auto get_seq_begin(const fqrec &rec) { return get_seq(rec); }
  [[nodiscard]] auto make_seq_begin(const fqrec &rec) { return get_seq(rec); }
  [[nodiscard]] auto get_seq_begin(const bamrec &) { return std::data(seq); }
  // clang-format on

  [[nodiscard]] auto
  make_seq_begin(const bamrec &rec) {
    const auto complement = [](const auto a) {
      return "TNGNNNCNNNNNNNNNNNNA"[a - 'A'];
    };
    auto rec_seq_itr = get_seq(rec);
    const auto rec_seq_end = get_seq_end(rec);
    if (rec.is_rev) {
      auto itr = std::begin(seq) + get_seq_size(rec);
      while (rec_seq_itr != rec_seq_end)
        *(--itr) = complement(*rec_seq_itr++);
    }
    else {
      auto itr = std::begin(seq);
      while (rec_seq_itr != rec_seq_end)
        *itr++ = *rec_seq_itr++;
    }
    return std::data(seq);
  }

  [[nodiscard]] auto
  process_quality_scores(const bamrec &rec) {
    return rec.is_rev
             ? count_quals_rev(get_qual(rec), get_qual_end(rec), qual_by_pos)
             : count_quals(get_qual(rec), get_qual_end(rec), qual_by_pos);
  }

  [[nodiscard]] auto
  process_quality_scores(const fqrec &rec) {
    return count_quals(get_qual(rec), get_qual_end(rec), qual_by_pos);
  }

  auto
  process_one_read_impl(const auto &rec) {
    // NOLINTBEGIN(cppcoreguidelines-pro-bounds-constant-array-index)
    static constexpr auto discrete_pct = [](const auto a, const auto b) {
      return (100 * a) / b;  // NOLINT(cppcoreguidelines-avoid-magic-numbers)
    };
    const auto read_len = static_cast<std::uint32_t>(get_seq_size(rec));
    if (read_len > max_read_len) {
      resize(read_len);
      if constexpr (std::is_same_v<std::decay_t<decltype(rec)>, bamrec>)
        seq.resize(read_len);
    }
    ++lengths[read_len];
    if (read_len == 0) [[unlikely]]
      return;
    max_read_len = read_len > max_read_len ? read_len : max_read_len;
    const auto seq_itr = make_seq_begin(rec);
    const auto seq_end = seq_itr + read_len;
    count_nucs(seq_itr, seq_end, base_counts);
    const auto gc = count_gc(seq_itr, seq_end);
    ++gc_content[discrete_pct(gc, read_len)];
    count_ns(seq_itr, seq_end, n_counts);
    const auto tot = process_quality_scores(rec);
    ++qual_by_read[tot / read_len];
    dr.count_seqs(seq_itr, read_len);
    am.match_adapters(seq_itr, read_len);
    // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  auto
  operator+=(const results_collector &rhs) -> const results_collector & {
    n_reads += rhs.n_reads;
    max_read_len = std::max(max_read_len, rhs.max_read_len);
    two_dim_add(base_counts, rhs.base_counts);
    add(gc_content, rhs.gc_content);
    vec_add(lengths, rhs.lengths);
    vec_add(n_counts, rhs.n_counts);
    two_dim_add(qual_by_pos, rhs.qual_by_pos);
    add(qual_by_read, rhs.qual_by_read);
    dr += rhs.dr;
    am += rhs.am;
    return *this;
  }

  template <typename self_t>
  [[nodiscard]] auto
  get_report(this self_t &self, const run_mode &mode, const file_info &info,
             const file_grades &grades) {
    const auto res = self.get_report_impl(mode, info, grades);
    std::string report;
    for (const auto &field : section_names)
      if (const auto res_itr = res.find(field); res_itr != std::cend(res))
        report += res_itr->second;
    return report;
  }

  [[nodiscard]] auto
  get_report_impl(const run_mode &mode, const file_info &info,
                  const file_grades &grades) const {
    const auto total_bases = tabular_dot(lengths);
    const auto gc_acc = [](const auto a, const auto &nuc) {
      return a + nuc[1] + nuc[3];  // NOLINT(*-avoid-magic-numbers)
    };
    const auto total_gc = std::accumulate(std::cbegin(base_counts),
                                          std::cend(base_counts), 0ul, gc_acc);
    const auto gt0 = [](const auto c) { return c > 0; };
    const std::uint64_t min_read_len =
      std::distance(std::cbegin(lengths), std::ranges::find_if(lengths, gt0));
    const auto groups = get_default_base_groups(max_read_len, do_groups(mode));
    return std::unordered_map<std::string, std::string>{
      {"basic_stats",
       format_basic_stats(info, n_reads, min_read_len, max_read_len, total_gc,
                          total_bases, grades)},
      {"quality_base", format_quality_base(qual_by_pos, groups, grades)},
      {"quality_sequence", format_quality_sequence(qual_by_read, grades)},
      {"sequence", format_sequence(base_counts, groups, grades)},
      {"gc_sequence", format_gc_sequence(gc_content, grades)},
      {"n_content", format_n_content(n_counts, base_counts, groups, grades)},
      {"sequence_length", format_sequence_length(lengths, grades)},
      {"duplication", dr.format_duplication(grades)},
      {"overrepresented", dr.format_overrepresented(grades)},
      {"adapter", am.get_report(n_reads, groups, grades)},
    };
  }

  template <typename self_t>
  [[nodiscard]] auto
  get_html(this self_t &self, const run_mode &mode, const file_info &info,
           const file_grades &grades) {
    const auto res = self.get_html_impl(mode, info, grades);
    std::string modules;
    for (const auto &section_name : section_names)
      if (const auto res_itr = res.find(section_name);
          res_itr != std::cend(res))
        modules += res_itr->second;
    return falco_get_html(info, grades, modules);
  }

  [[nodiscard]] auto
  get_html_impl(const run_mode &mode, const file_info &info,
                const file_grades &grades) const {
    const auto total_bases = tabular_dot(lengths);
    const auto gc_acc = [](const auto a, const auto &nuc) {
      return a + nuc[1] + nuc[3];  // NOLINT(*-avoid-magic-numbers)
    };
    const auto total_gc = std::accumulate(std::cbegin(base_counts),
                                          std::cend(base_counts), 0ul, gc_acc);
    const auto gt0 = [](const auto c) { return c > 0; };
    const std::uint64_t min_read_len =
      std::distance(std::cbegin(lengths), std::ranges::find_if(lengths, gt0));
    const auto groups = get_default_base_groups(max_read_len, do_groups(mode));
    const auto basic_stats_html = format_basic_stats_html(
      info, n_reads, min_read_len, max_read_len, total_gc, total_bases, grades);
    return std::unordered_map<std::string, std::string>{
      {"basic_stats", basic_stats_html},
      {"quality_base", format_quality_base_html(qual_by_pos, groups, grades)},
      {"quality_sequence", format_quality_sequence_html(qual_by_read, grades)},
      {"sequence", format_sequence_html(base_counts, groups, grades)},
      {"gc_sequence", format_gc_sequence_html(gc_content, grades)},
      {"n_content",
       format_n_content_html(n_counts, base_counts, groups, grades)},
      {"sequence_length", format_sequence_length_html(lengths, grades)},
      {"duplication", dr.format_duplication_html(grades)},
      {"overrepresented", dr.format_overrepresented_html(grades)},
      {"adapter", am.get_html(n_reads, groups, grades)},
    };
  }
};

struct results_collector_tile : public results_collector {
  tile_processor tp;
  bool has_tiles{};

  auto
  init_impl(const auto &info) {
    results_collector::init_impl(info);
    has_tiles = info.has_tiles;
    if (has_tiles)
      tp.init(info);
  };

  auto
  get_grades_impl(file_grades &grades) {
    results_collector::get_grades_impl(grades);
    if (has_tiles)
      grades.emplace("tile", tp.get_grade());
  }

  auto
  finalize_impl(const run_mode &mode, const file_info &info) {
    results_collector::finalize_impl(mode, info);
    if (has_tiles)
      tp.finalize(mode, info);
  }

  auto
  process_one_read_impl(const auto &rec) {
    results_collector::process_one_read_impl(rec);
    if (has_tiles) [[likely]]
      tp(rec);
  }

  auto
  operator+=(const results_collector_tile &rhs)
    -> const results_collector_tile & {
    results_collector::operator+=(rhs);
    tp += rhs.tp;
    return *this;
  }

  [[nodiscard]] auto
  get_report_impl(const run_mode &mode, const file_info &info,
                  const file_grades &grades) const {
    assert(array_contains(section_names, "tile"));
    auto report = results_collector::get_report_impl(mode, info, grades);
    if (has_tiles) {
      const auto groups =
        get_default_base_groups(max_read_len, do_groups(mode));
      report.emplace("tile", tp.get_report(groups, grades.grade("tile")));
    }
    return report;
  }

  [[nodiscard]] auto
  get_html_impl(const run_mode &mode, const file_info &info,
                const file_grades &grades) const {
    assert(array_contains(section_names, "tile"));
    auto html = results_collector::get_html_impl(mode, info, grades);
    if (has_tiles) {
      const auto groups =
        get_default_base_groups(max_read_len, do_groups(mode));
      html.emplace("tile", tp.get_html(groups, grades));
    }
    return html;
  }
};

struct results_collector_kmer : public results_collector {
  kmer_counter kc;

  auto
  get_grades_impl(file_grades &grades) {
    results_collector::get_grades_impl(grades);
    grades.emplace("kmer", kc.get_grade());
  }

  auto
  finalize_impl(const run_mode &mode, const file_info &info) {
    results_collector::finalize_impl(mode, info);
    kc.finalize(mode);
  }

  auto
  process_one_read_impl(const auto &rec) {
    results_collector::process_one_read_impl(rec);
    kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
  }

  auto
  operator+=(const results_collector_kmer &rhs)
    -> const results_collector_kmer & {
    results_collector::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  get_report_impl(const run_mode &mode, const file_info &info,
                  const file_grades &grades) const {
    assert(array_contains(section_names, "kmer"));
    auto report = results_collector::get_report_impl(mode, info, grades);
    report.emplace("kmer", kc.get_report(grades.grade("kmer")));
    return report;
  }

  [[nodiscard]] auto
  get_html_impl(const run_mode &mode, const file_info &info,
                const file_grades &grades) const {
    assert(array_contains(section_names, "kmer"));
    auto html = results_collector::get_html_impl(mode, info, grades);
    html.emplace("kmer", kc.get_html(grades));
    return html;
  }
};

struct results_collector_tile_kmer : public results_collector_tile {
  kmer_counter kc;

  auto
  get_grades_impl(file_grades &grades) {
    results_collector_tile::get_grades_impl(grades);
    grades.emplace("kmer", kc.get_grade());
  }

  auto
  finalize_impl(const run_mode &mode, const file_info &info) {
    results_collector_tile::finalize_impl(mode, info);
    kc.finalize(mode);
  }

  auto
  process_one_read_impl(const auto &rec) {
    results_collector_tile::process_one_read_impl(rec);
    kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
  }

  auto
  operator+=(const results_collector_tile_kmer &rhs)
    -> const results_collector_tile_kmer & {
    results_collector_tile::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  get_report_impl(const run_mode &mode, const file_info &info,
                  const file_grades &grades) const {
    assert(array_contains(section_names, "kmer"));
    auto report = results_collector_tile::get_report_impl(mode, info, grades);
    report.emplace("kmer", kc.get_report(grades.grade("kmer")));
    return report;
  }

  [[nodiscard]] auto
  get_html_impl(const run_mode &mode, const file_info &info,
                const file_grades &grades) const {
    assert(array_contains(section_names, "kmer"));
    auto html = results_collector_tile::get_html_impl(mode, info, grades);
    html.emplace("kmer", kc.get_html(grades));
    return html;
  }
};

template <typename results_t>
static inline auto
accumulate_results(std::vector<results_t> &r, const auto file_id) {
  // pointer jumping strategy but only for the specified file
  auto id = std::views::iota(0, std::ssize(r)) | std::ranges::to<std::vector>();
  while (std::size(id) > 1) {
    {
      std::vector<std::jthread> workers;
      for (auto i = 0u; i + 1 < std::size(id); i += 2)
        // NOLINTNEXTLINE(performance-inefficient-vector-operation)
        workers.emplace_back(
          [&, i] { r[id[i]][file_id] += r[id[i + 1]][file_id]; });
    }
    auto j = 0;
    for (auto i = 0; i + 1 < std::ssize(id); i += 2)
      id[j++] = id[i];
    if (std::ssize(id) % 2 == 1)
      id[j++] = id.back();
    id.resize(j);
  }
}

#endif  // SRC_RESULTS_COLLECTOR_HPP_
