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

#ifndef SRC_FALCO_RESULTS_HPP_
#define SRC_FALCO_RESULTS_HPP_

#include "adapter_matcher.hpp"
#include "bam_file.hpp"
#include "duplication_results.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "format_output.hpp"
#include "kmer_counter.hpp"
#include "tile_processor.hpp"

#include <algorithm>
#include <cstdint>
#include <ranges>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

// clang-format off
static constexpr auto report_section_order = std::array{
  "basic_stats",
  "qual_by_pos",
  "qual_by_read",
  "base_composition",
  "gc_content",
  "n_counts",
  "read_lengths",
  "duplication_levels",
  "overrepresented",
  "adapters",
  "tiles",
  "kmers",
};
// clang-format on

[[nodiscard]] static inline auto
array_contains(const auto &a, const auto &b) {
  return std::ranges::find(a, b) != std::cend(a);
};

static constexpr auto assumed_page_size = 4096;
struct alignas(assumed_page_size) falco_results {
  // struct alignas(std::hardware_destructive_interference_size) falco_results {
  std::uint64_t n_reads{};
  std::uint64_t max_read_len{};
  std::vector<falco::nuc_array> nucs;
  falco::gc_content_array gcs{};
  std::vector<std::uint64_t> n_counts;
  std::vector<std::uint64_t> lengths;
  std::vector<falco::qual_array> qual_by_pos;
  falco::qual_array qual_by_read{};
  duplication_results dr;
  adapter_matcher am;
  std::string seq;

  falco_results() : lengths(1, 0) {}  // in case all reads have length 0

  auto
  resize(const std::uint32_t updated_length) {
    nucs.resize(updated_length);
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
  finalize(this self_t &self, const run_mode &mode, const file_info &info) {
    self.finalize_impl(mode, info);
  }

  auto
  finalize_impl(const run_mode &mode, const file_info &info) {
    if (!is_mapped_reads(info.format))
      adjust_fastq_qual_encoding(qual_by_pos, qual_by_read, info.encoding);
    am.finalize(mode);
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
    count_nucs(seq_itr, seq_end, nucs);
    const auto gc = count_gc(seq_itr, seq_end);
    ++gcs[discrete_pct(gc, read_len)];
    count_ns(seq_itr, seq_end, n_counts);
    const auto tot = process_quality_scores(rec);
    ++qual_by_read[tot / read_len];
    dr.count_seqs(seq_itr, read_len);
    am.match_adapters(seq_itr, read_len);
    // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  auto
  operator+=(const falco_results &rhs) -> const falco_results & {
    n_reads += rhs.n_reads;
    max_read_len = std::max(max_read_len, rhs.max_read_len);
    two_dim_add(nucs, rhs.nucs);
    add(gcs, rhs.gcs);
    vec_add(lengths, rhs.lengths);
    vec_add(n_counts, rhs.n_counts);
    two_dim_add(qual_by_pos, rhs.qual_by_pos);
    add(qual_by_read, rhs.qual_by_read);
    dr += rhs.dr;
    am += rhs.am;
    return *this;
  }

  [[nodiscard]] auto
  adjust_nucs_for_ns() const {  // Ns were counted among the C in nucs
    auto nucs_no_n = nucs;
    for (auto i = 0u; i < std::size(nucs_no_n); ++i)
      nucs_no_n[i][3] -= n_counts[i];
    return nucs_no_n;
  };

  template <typename self_t>
  [[nodiscard]] auto
  get_report(this self_t &self, const file_info &info, grades &g) {
    const auto res = self.get_report_impl(info, g);
    std::string report;
    for (const auto &field : report_section_order)
      if (const auto res_itr = res.find(field); res_itr != std::cend(res))
        report += res_itr->second;
    return report;
  }

  [[nodiscard]] auto
  get_report_impl(const file_info &info, grades &g) const {
    const auto nucs_no_n = adjust_nucs_for_ns();
    const auto total_nucs = tabular_dot(lengths);
    const auto gc_acc = [](const auto a, const auto &nuc) {
      return a + nuc[1] + nuc[3];  // NOLINT(*-avoid-magic-numbers)
    };
    const auto total_gc = std::accumulate(std::cbegin(nucs_no_n),
                                          std::cend(nucs_no_n), 0ul, gc_acc);
    const auto gt0 = [](const auto c) { return c > 0; };
    const std::uint64_t min_read_len =
      std::distance(std::cbegin(lengths), std::ranges::find_if(lengths, gt0));
    const auto encoding_label = get_quality_score_label(info.encoding);
    return std::map<std::string, std::string>{
      {"basic_stats",
       format_basic_stats(info.name, n_reads, min_read_len, max_read_len,
                          total_gc, total_nucs, encoding_label, g.basic_stats)},
      {"qual_by_pos", format_qual_by_pos(qual_by_pos, g.qual_by_pos)},
      {"qual_by_read", format_qual_by_read(qual_by_read, g.qual_by_read)},
      {"base_composition",
       format_base_composition(nucs_no_n, g.base_composition)},
      {"gc_content", format_gc_content(gcs, g.gc_content)},
      {"n_counts", format_n_counts(n_counts, nucs, g.n_counts)},
      {"read_lengths", format_read_lengths(lengths, g.read_lengths)},
      {"duplication_levels",
       dr.format_duplication_levels(g.duplication_levels)},
      {"overrepresented", dr.format_overrepresented(g.overrepresented)},
      {"adapters", am.get_report(n_reads, g.adapter_content)},
    };
    // auto r =
    //   format_basic_stats(info.name, n_reads, min_read_len, max_read_len,
    //                      total_gc, total_nucs, encoding_label,
    //                      g.basic_stats);
    // r += format_qual_by_pos(qual_by_pos, g.qual_by_pos);
    // r += format_qual_by_read(qual_by_read, g.qual_by_read);
    // r += format_base_composition(nucs_no_n, g.base_composition);
    // r += format_gc_content(gcs, g.gc_content);
    // r += format_n_counts(n_counts, nucs, g.n_counts);
    // r += format_read_lengths(lengths, g.read_lengths);
    // r += dr.format_duplication_levels(g.duplication_levels);
    // r += dr.format_overrepresented(g.overrepresented);
    // r += am.get_report(n_reads, g.adapter_content);
    // return r;

    // clang-format on
  }

  template <typename self_t>
  [[nodiscard]] auto
  get_html(this self_t &self, const file_info &info) {
    return self.get_html_impl(info);
  }

  [[nodiscard]] auto
  get_html_impl([[maybe_unused]] const file_info &info) const {
#ifdef MAKE_HTML
    const auto nucs_no_n = adjust_nucs_for_ns();
    const auto total_nucs = tabular_dot(lengths);
    const auto gc_acc = [](const auto a, const auto &nuc) {
      return a + nuc[1] + nuc[3];  // NOLINT(*-avoid-magic-numbers)
    };
    const auto total_gc = std::accumulate(std::cbegin(nucs_no_n),
                                          std::cend(nucs_no_n), 0ul, gc_acc);
    const auto gt0 = [](const auto c) { return c > 0; };
    const std::uint64_t min_read_len =
      std::distance(std::cbegin(lengths), std::ranges::find_if(lengths, gt0));
    const auto encoding_label = get_quality_score_label(info.encoding);
    auto r =
      format_basic_stats_html(info.name, n_reads, min_read_len, max_read_len,
                              total_gc, total_nucs, encoding_label);
    r += format_qual_by_pos_html(qual_by_pos);
    r += format_qual_by_read_html(qual_by_read);
    r += format_base_composition_html(nucs);
    r += format_n_counts_html(n_counts, nucs);
    r += format_gc_content_html(gcs);
    r += format_read_lengths_html(lengths);
    r += dr.format_duplication_levels_html();
    r += dr.format_overrepresented_html();
    return r;
#else
    return std::string{"not yet implemented: should be ready soon!"};
#endif  // MAKE_HTML
  }
};

struct falco_results_tile : public falco_results {
  tile_processor tp;
  bool has_tiles{};

  auto
  init_impl(const auto &info) {
    falco_results::init_impl(info);
    has_tiles = info.has_tiles;
    tp.init(info);
  };

  auto
  finalize_impl(const run_mode &mode, const file_info &info) {
    falco_results::finalize_impl(mode, info);
    tp.finalize(mode, info);
  }

  auto
  process_one_read_impl(const auto &rec) {
    falco_results::process_one_read_impl(rec);
    if (has_tiles) [[likely]]
      tp(rec);
  }

  auto
  operator+=(const falco_results_tile &rhs) -> const falco_results_tile & {
    falco_results::operator+=(rhs);
    tp += rhs.tp;
    return *this;
  }

  [[nodiscard]] auto
  get_report_impl(const file_info &info, grades &g) const {
    assert(array_contains(report_section_order, "tiles"));
    auto report = falco_results::get_report_impl(info, g);
    report.emplace("tiles", tp.get_report(g.tile_analaysis));
    return report;
  }

  [[nodiscard]] auto
  get_html_impl(const file_info &info) const {
    return falco_results::get_html_impl(info) + tp.get_html();
  }
};

struct falco_results_kmer : public falco_results {
  kmer_counter kc;

  auto
  process_one_read_impl(const auto &rec) {
    falco_results::process_one_read_impl(rec);
    kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
  }

  auto
  operator+=(const falco_results_kmer &rhs) -> const falco_results_kmer & {
    falco_results::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  get_report_impl(const file_info &info, grades &g) const {
    assert(array_contains(report_section_order, "kmers"));
    auto report = falco_results::get_report_impl(info, g);
    report.emplace("kmers", kc.get_report(g.kmer_content));
    return report;
  }

  [[nodiscard]] auto
  get_html_impl(const file_info &info) const {
    return falco_results::get_html_impl(info) + kc.get_html();
  }
};

struct falco_results_tile_kmer : public falco_results_tile {
  kmer_counter kc;

  auto
  finalize_impl(const run_mode &mode, const file_info &info) {
    falco_results_tile::finalize_impl(mode, info);
  }

  auto
  process_one_read_impl(const auto &rec) {
    falco_results_tile::process_one_read_impl(rec);
    kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
  }

  auto
  operator+=(const falco_results_tile_kmer &rhs)
    -> const falco_results_tile_kmer & {
    falco_results_tile::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  get_report_impl(const file_info &info, grades &g) const {
    assert(array_contains(report_section_order, "kmers"));
    auto report = falco_results_tile::get_report_impl(info, g);
    report.emplace("kmers", kc.get_report(g.kmer_content));
    return report;
  }

  [[nodiscard]] auto
  get_html_impl(const file_info &info) const {
    return falco_results_tile::get_html_impl(info) + kc.get_html();
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

#endif  // SRC_FALCO_RESULTS_HPP_
