// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_DUPLICATION_RESULTS_HPP_
#define SRC_DUPLICATION_RESULTS_HPP_

#include "falco_word.hpp"

#include "boost/boost_unordered.hpp"

#include <cstdint>
#include <iterator>
#include <numeric>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

class run_mode;
struct file_grades;
struct file_info;

using dups_map_t = boost::unordered_flat_map<falco_word, std::uint64_t>;

struct overrep_t {
  falco_word w;
  std::uint64_t n_obs{};     // number of observations
  double pct_val{};          // corresponding percentage
  std::int32_t contam_id{};  // index of matching contaminant
};

struct dup_summary_t {
  std::uint64_t max_dup{};
  std::uint64_t n_reads{};
  std::vector<std::uint64_t> hist_mass;
  std::vector<std::uint64_t> hist_dedup;
};

struct dups_init_t {
  std::int64_t count_at_limit{};
  dups_map_t dups_zero;
  dups_init_t() = default;
  dups_init_t(const dups_map_t &dups) {
    const auto vals = dups | std::views::values;
    count_at_limit = static_cast<std::int64_t>(
      std::reduce(std::cbegin(vals), std::cend(vals)));
    for (const auto &fw : dups | std::views::keys)
      dups_zero.emplace(fw, 0);
  }
  operator bool() const { return !dups_zero.empty(); }
};

struct duplication_results {
  static constexpr auto max_n_reads_total{1'000'000};
  static constexpr auto max_reads_to_hash{100'000};
  static constexpr auto default_read_skip{10};
  static constexpr auto overrep_cutoff = 0.001;

  std::int64_t count_at_limit{};
  std::int64_t read_skip{default_read_skip};
  std::int64_t read_idx{};
  dups_map_t dups;

  auto
  initialize(const run_mode &mode, const file_info &info,
             const dups_init_t &dups_init) -> void;

  auto
  initialize(const run_mode &mode, const file_info &info) -> void;

  auto
  get_n_counted_reads() const -> std::uint64_t;

  auto
  release() -> void {
    dups.clear();
    dups_map_t().swap(dups);
  }

  [[nodiscard]] auto
  get_preseq_hist() const
    -> std::vector<std::pair<std::uint32_t, std::uint32_t>>;

  [[nodiscard]] auto
  get_dups_summary(const std::uint64_t n_reads) const -> dup_summary_t;

  [[nodiscard]] auto
  get_dups_summary() const -> dup_summary_t;

  [[nodiscard]] auto
  get_overrepresented(const std::uint64_t n_reads) const
    -> std::vector<overrep_t>;

  auto
  add_and_consume(duplication_results &rhs) -> void;

  auto
  count_seqs_orig(const auto seq_itr, const auto sz) {
    auto itr = dups.find(falco_word(seq_itr, sz));
    if (itr != std::cend(dups))
      ++itr->second;
  }

  auto
  count_seqs(const auto seq_itr, const auto sz) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      ++dups[falco_word(seq_itr, sz)];
    }
  }
};

[[nodiscard]] auto
get_grade_duplication(const dup_summary_t &summary) -> std::string;

[[nodiscard]] auto
duplication_report(const dup_summary_t &summary,
                   const file_grades &grades) -> std::string;

[[nodiscard]] auto
duplication_html(const dup_summary_t &summary,
                 const file_grades &grades) -> std::string;

[[nodiscard]] auto
get_grade_overrepresented(const std::uint64_t n_reads,
                          const duplication_results &dr) -> std::string;

[[nodiscard]] auto
overrepresented_report(const std::vector<overrep_t> &overrep,
                       const file_grades &grades) -> std::string;

[[nodiscard]] auto
overrepresented_html(const std::vector<overrep_t> &overrep,
                     const file_grades &grades) -> std::string;

#endif  // SRC_DUPLICATION_RESULTS_HPP_
