// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "original_duplicates.hpp"

#include "duplication_results.hpp"
#include "falco_file_format.hpp"
#include "file_info.hpp"

#include "boost/boost_unordered.hpp"

#include <htslib/bgzf.h>  // IWYU pragma: keep
#include <htslib/hfile.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <mutex>
#include <ranges>
#include <string>
#include <thread>
#include <utility>  // for std::move
#include <vector>

[[nodiscard]] auto
init_orig_dups_bam(const std::string &filename, const std::uint64_t n_unique)
  -> dups_map_t {
  // ADS: n_unique_multiplier: contributes to max reads to examine, multiplying
  // the n_unique that is the target; impossible to know how many needed to hit
  // the target, so try to overshoot.
  static constexpr auto n_unique_multiplier = 10;
  static constexpr auto complem = [](char &x) {
    assert(x >= 'A' && x <= 'T');
    x = "TNGNNNCNNNNNNNNNNNNA"[x - 'A'];
  };
  static constexpr auto to_acgtn = [](const char x) {
    static constexpr char mask = 15;
    static constexpr std::span seq_nt16_str_acgtn = "NACNGNNNTNNNNNNN";
    return seq_nt16_str_acgtn[x & mask];
  };
  const std::unique_ptr<htsFile, int (*)(htsFile *)> in(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!in)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
  const std::unique_ptr<sam_hdr_t, void (*)(sam_hdr_t *)> h(
    sam_hdr_read(in.get()), &sam_hdr_destroy);
  if (!h)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to read header: " + filename);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-union-access)
  std::unique_ptr<bam1_t, void (*)(bam1_t *)> rec(bam_init1(), &bam_destroy1);

  dups_map_t dups;
  std::vector<char> buffer;
  int r{};
  std::uint64_t n_reads{};
  const auto max_n_reads = n_unique_multiplier * n_unique;
  while (n_reads++ < max_n_reads && std::size(dups) < n_unique &&
         (r = sam_read1(in.get(), h.get(), rec.get())) >= 0) {
    const auto l_qseq = rec->core.l_qseq;
    const auto seq = bam_get_seq(rec.get());
    if (std::ssize(buffer) < l_qseq)
      buffer.resize(l_qseq);
    for (auto i = 0; i < l_qseq; ++i)
      buffer[i] = to_acgtn(bam_seqi(seq, i));
    if (bam_is_rev(rec.get())) {
      std::reverse(std::begin(buffer), std::begin(buffer) + l_qseq);
      std::transform(std::cbegin(buffer), std::cbegin(buffer) + l_qseq,
                     std::begin(buffer), complem);
    }
    ++dups[falco_word(std::data(buffer), l_qseq)];
  }
  if (r < -1)  // error
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "error reading bam record from: " + filename);
  return dups;
}

[[nodiscard]] auto
init_orig_dups_fq(const std::string &filename, const std::uint64_t n_unique)
  -> dups_map_t {
  // ADS: n_unique_multiplier: contributes to max reads to examine, multiplying
  // the n_unique that is the target; impossible to know how many needed to hit
  // the target, so try to overshoot.
  static constexpr auto n_unique_multiplier = 10;
  static constexpr auto to_acgtn = [](const auto c) {
    static constexpr auto mask = 31;
    static constexpr std::span to_acgtn_sp = "NANCNNNGNNNNNNNNNNNNTNNNNNNNNNNN";
    return to_acgtn_sp[c & mask];
  };
  const std::unique_ptr<BGZF, int (*)(BGZF *)> in(
    bgzf_open(std::data(filename), "r"), &bgzf_close);
  if (!in)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
  dups_map_t dups;
  kstring_t s = KS_INITIALIZE;
  int r{};
  std::uint64_t n_lines{};
  std::uint64_t n_reads{};
  const auto max_n_reads = n_unique_multiplier * n_unique;
  while (n_reads < max_n_reads && std::size(dups) < n_unique &&
         (r = bgzf_getline(in.get(), '\n', &s)) >= 0) {
    if (n_lines % 4 == 1) {
      const auto l_qseq = ks_len(&s);
      std::ranges::transform(std::span(s.s, l_qseq), s.s, to_acgtn);
      const auto seq = ks_c_str(&s);
      ++dups[falco_word(seq, l_qseq)];
      ++n_reads;
    }
    ++n_lines;
  }
  if (r < -1)  // error
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "error reading bam record from: " + filename);
  ks_free(&s);
  return dups;
}

[[nodiscard]] auto
initialize_original_duplicates(const std::vector<std::string> &infiles,
                               const std::vector<file_info> &infos,
                               const std::uint32_t n_threads)
  -> std::vector<dups_init_t> {
  assert(std::size(infiles) == std::size(infos));
  const std::uint64_t n_files = std::size(infiles);
  std::vector<dups_map_t> dups(n_files);
  auto n_active_files = n_files;
  std::vector<std::jthread> workers;
  const auto n_workers =
    std::min(n_files, static_cast<std::uint64_t>(n_threads));
  workers.reserve(n_workers);
  std::mutex mtx;
  // NOLINTNEXTLINE(clang-analyzer-deadcode.DeadStores)
  for (const auto _ : std::views::iota(0u, n_workers))
    workers.emplace_back([&] {
      std::uint64_t file_id{};
      while (true) {
        {
          const std::scoped_lock l(mtx);
          if (n_active_files == 0)
            return;
          file_id = --n_active_files;
        }
        dups[file_id] =
          falco::is_mapped_reads(infos[file_id].format)
            ? init_orig_dups_bam(infiles[file_id],
                                 duplication_results::max_reads_to_hash)
            : init_orig_dups_fq(infiles[file_id],
                                duplication_results::max_reads_to_hash);
      }
    });
  std::ranges::for_each(workers, [](auto &w) { w.join(); });
  std::vector<dups_init_t> ret;
  ret.reserve(n_files);
  std::ranges::for_each(dups,
                        [&](auto &&d) { ret.emplace_back(std::move(d)); });
  return ret;
}
