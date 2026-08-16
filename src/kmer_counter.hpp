// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_KMER_COUNTER_HPP_
#define SRC_KMER_COUNTER_HPP_

#include "falco_utils.hpp"

#include <array>
#include <cstdint>
#include <iterator>
#include <string>
#include <thread>  // IWYU pragma: keep
#include <vector>

class run_mode;

struct kmer_result {
  // holds info for {kmer} x {position} to be sorted, filtered and output
  std::uint64_t kmer{};
  std::uint64_t count{};
  double pval{1.0};
  double obs_exp{};
  std::uint64_t pos{};
  [[nodiscard]] auto
  string() const -> std::string;
  [[nodiscard]] auto
  decode() const -> std::string;
};

struct kmer_counter {
  // cutoffs for what to report
  static constexpr auto min_obs_exp_to_report = 5.0;
  static constexpr auto max_pval_to_report = 0.01;
  [[maybe_unused]] static constexpr auto max_kmers_to_plot = 10;
  static constexpr auto n_kmers_to_report = 20;

  // don't count kmers past this length of read prefix
  static constexpr auto max_len_for_analysis = 500;

  static constexpr auto read_skip = 50 - 1;
  static constexpr auto kmer_size = 7;
  static constexpr auto n_kmers = ipow(4, kmer_size);
  static constexpr auto kmer_mask = n_kmers - 1;

  std::int32_t read_idx{};
  std::int32_t max_read_len{};
  std::vector<std::array<std::uint64_t, n_kmers>> kmer_counts;

  auto
  resize(const std::int32_t updated_length) {
    max_read_len = updated_length;
    kmer_counts.resize(updated_length);
  }

  auto
  release() -> void {
    kmer_counts.clear();
    kmer_counts.shrink_to_fit();
  }

  static constexpr auto
  shift_kmer(auto &k, auto c) {
    k = (k << 2) | encode(c);
  };

  auto
  count_kmers(auto seq_itr, auto sz) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      sz = sz > max_len_for_analysis ? max_len_for_analysis : sz;
      if (sz > static_cast<decltype(sz)>(max_read_len)) [[unlikely]]
        resize(sz);
      const auto seq_end = seq_itr + sz;
      auto out_itr = std::begin(kmer_counts);
      auto kmer = 0;
      for (auto i = 0; i < kmer_size - 1 && seq_itr != seq_end; ++i)
        shift_kmer(kmer, *seq_itr++);
      while (seq_itr != seq_end) {
        shift_kmer(kmer, *seq_itr++);
        // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-constant-array-index)
        ++(*out_itr++)[kmer & kmer_mask];
      }
    }
  }

  auto
  apply_groups(const run_mode &mode) -> void;

  auto
  operator+=(const kmer_counter &rhs) -> const kmer_counter &;

  auto
  add_and_consume(kmer_counter &&rhs) -> void;

  [[nodiscard]] auto
  get_kmer_results() const -> std::vector<kmer_result>;

  [[nodiscard]] static auto
  decode_kmer(auto word, const auto n_bases) -> std::string;
};

[[nodiscard]] auto
get_grade_kmer(const std::vector<kmer_result> &results) -> std::string;

#endif  // SRC_KMER_COUNTER_HPP_
