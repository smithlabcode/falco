// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_RESULTS_COLLECTOR_HPP_
#define SRC_RESULTS_COLLECTOR_HPP_

#include "adapter_matcher.hpp"
#include "bamrec.hpp"
#include "duplication_results.hpp"
#include "falco_file_format.hpp"
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "file_info.hpp"
#include "fqrec.hpp"
#include "kmer_counter.hpp"
#include "quality_score.hpp"
#include "run_mode.hpp"
#include "samrec.hpp"
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
}

static constexpr auto assumed_page_size = 4096;
struct alignas(assumed_page_size) results_collector {
  // alignas(std::hardware_destructive_interference_size)
  std::uint64_t n_reads{};
  std::uint64_t max_read_len{};
  std::vector<falco::nuc_array> base_counts;
  std::vector<falco::gc_content_array> gc_content;
  std::vector<std::uint64_t> n_counts;
  std::vector<std::uint64_t> lengths;
  std::vector<falco::qual_array> qual_by_pos;
  falco::qual_array qual_by_read{};
  duplication_results dr;
  adapter_matcher am;
  tile_processor tp;
  kmer_counter kc;
  bool do_tiles{};
  bool do_kmers{};

  results_collector() : lengths(1, 0) {}  // in case all reads have length 0

  auto
  resize(const std::uint32_t updated_length) {
    base_counts.resize(updated_length);
    if (std::size(gc_content) < falco::gc_content_array_max_size)
      resize_gc_content(updated_length, gc_content);
    n_counts.resize(updated_length);
    lengths.resize(updated_length + 1);  // need one extra here
    qual_by_pos.resize(updated_length);
    am.resize(updated_length);
  }

  auto
  init(const run_mode &mode, const auto &info) {
    dr.initialize(mode, info);
    do_tiles = mode.do_tiles() && info.has_tiles;
    do_kmers = mode.do_kmers();
    if (do_tiles)
      tp.init(info);
  }

  auto
  adjust_base_counts_for_ns() -> void {
    // subtract N counts from the position where they were counted
    for (auto i = 0u; i < std::size(base_counts); ++i)
      base_counts[i][unknown_base_index] -= n_counts[i];
  }

  auto
  process_reads_bam(bamrec::pos_t cursor, const bamrec::pos_t lim) -> void {
    bamrec rec{};
    while (cursor < lim && get_next(cursor, lim, rec)) {
      process_one_read(rec);
      ++n_reads;
    }
  }

  auto
  process_reads_sam(samrec::pos_t cursor, const samrec::pos_t lim) -> void {
    samrec rec{};
    while (cursor < lim && get_next(cursor, lim, rec)) {
      process_one_read(rec);
      ++n_reads;
    }
  }

  auto
  process_reads_fq(fqrec::pos_t cursor, const fqrec::pos_t lim) -> void {
    fqrec rec{};
    while (cursor < lim && (rec = get_next(cursor, lim))) {
      process_one_read(rec);
      ++n_reads;
    }
  }

  auto
  process_one_read(const auto &rec) -> void {
    // NOLINTBEGIN(cppcoreguidelines-pro-bounds-constant-array-index)
    static constexpr auto get_arr_idx = [&](const std::integral auto l) {
      constexpr auto x =
        static_cast<decltype(l)>(falco::gc_content_array_max_lim);
      return l < x ? l : x;
    };
    static constexpr auto get_gc_idx = [&](const std::integral auto a,
                                           const std::integral auto b) {
      constexpr auto x =
        static_cast<decltype(b)>(falco::gc_content_array_max_lim);
      return b < x ? a : (x * a) / b;
    };
    const auto read_len = static_cast<std::uint32_t>(get_seq_size(rec));
    if (read_len > max_read_len) [[unlikely]]
      resize(read_len);
    ++lengths[read_len];
    if (read_len == 0) [[unlikely]]
      return;
    max_read_len = read_len > max_read_len ? read_len : max_read_len;
    const auto seq_itr = get_seq(rec);
    const auto seq_end = get_seq_end(rec);
    count_nucs(seq_itr, seq_end, base_counts);
    const auto gc = count_gc(seq_itr, seq_end);
    assert(get_arr_idx(read_len) < std::ssize(gc_content));
    ++gc_content[get_arr_idx(read_len)][get_gc_idx(gc, read_len)];
    count_ns(seq_itr, seq_end, n_counts);
    const auto tot = count_quals(get_qual(rec), get_qual_end(rec), qual_by_pos);
    ++qual_by_read[tot / read_len];
    dr.count_seqs(seq_itr, read_len);
    am.match_adapters(seq_itr, read_len);
    if (do_tiles)
      tp(rec);
    if (do_kmers)
      kc.count_kmers(seq_itr, read_len);
    // NOLINTEND(cppcoreguidelines-pro-bounds-constant-array-index)
  }

  auto
  add_and_consume(results_collector &&rhs) -> void {
    n_reads += rhs.n_reads;
    max_read_len = std::max(max_read_len, rhs.max_read_len);
    two_dim_add_and_consume(base_counts, rhs.base_counts);
    two_dim_add_and_consume(gc_content, std::move(rhs.gc_content));
    vec_add_and_consume(lengths, std::move(rhs.lengths));
    vec_add_and_consume(n_counts, std::move(rhs.n_counts));
    two_dim_add(qual_by_pos, rhs.qual_by_pos);
    add(qual_by_read, rhs.qual_by_read);
    dr.add_and_consume(rhs.dr);
    am.add_and_consume(std::move(rhs.am));
    if (do_tiles)
      tp.add_and_consume(std::move(rhs.tp));
    if (do_kmers)
      kc.add_and_consume(std::move(rhs.kc));
  }

  auto
  finalize(file_info &info) {
    adjust_base_counts_for_ns();
    info.set_encoding(identify_encoding(qual_by_pos, info));
    if (!is_bam(info.format))
      adjust_fastq_qual_encoding(qual_by_pos, qual_by_read, info.encoding);
    if (do_tiles)
      tp.finalize(info);
  }
};

inline auto
accumulate_results(std::vector<std::vector<results_collector>> &r,
                   const std::int32_t file_id) {
  auto j = 1u;
  while (j < std::size(r)) {
    std::vector<std::jthread> workers;
    workers.reserve(std::size(r) / (2 * j));
    for (auto i = 0u; i + j < std::size(r); i += 2 * j)
      workers.emplace_back([&, i, j] {
        r[i][file_id].add_and_consume(std::move(r[i + j][file_id]));
      });
    j *= 2;
  }
}

inline auto
process_reads(results_collector &results, fq_task_t &task) {
  results.process_reads_fq(task.beg, task.end);
}

inline auto
process_reads(results_collector &results, bam_task_t &task) {
  results.process_reads_bam(task.beg, task.end);
}

inline auto
process_reads(results_collector &results, sam_task_t &task) {
  results.process_reads_sam(task.beg, task.end);
}

#endif  // SRC_RESULTS_COLLECTOR_HPP_
