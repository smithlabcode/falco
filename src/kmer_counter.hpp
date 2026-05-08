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

#ifndef SRC_KMER_COUNTER_HPP_
#define SRC_KMER_COUNTER_HPP_

#include "falco_utils.hpp"

#include <array>
#include <cstdint>
#include <format>
#include <iterator>
#include <string>
#include <thread>  // IWYU pragma: keep
#include <utility>
#include <vector>

struct kmer_counter {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.02, "pass"},
    std::pair{0.05, "warn"},
    std::pair{1.00, "error"},
  };

  static constexpr auto kmer_step = 50;
  static constexpr auto kmer_size = 7;
  static constexpr auto n_kmers = ipow(4, kmer_size);
  static constexpr auto kmer_mask = n_kmers - 1;
  static constexpr auto max_kmers_to_report = 20;
  static constexpr auto max_kmers_to_plot = 10;

  std::uint64_t next_kmer_read{};
  std::vector<std::array<std::uint64_t, n_kmers>> kmer_counts;

  explicit kmer_counter(const std::uint32_t init_read_len) :
    kmer_counts(init_read_len, std::array<std::uint64_t, n_kmers>{}) {}

  auto
  resize(const std::uint32_t updated_length) {
    kmer_counts.resize(updated_length);
  }

  static constexpr auto
  shift_kmer(auto &k, auto c) {
    k = ((k << 2) | encode(c)) & kmer_mask;
  };

  auto
  count_kmers(auto seq_itr, const auto sz) {
    const auto seq_end = seq_itr + sz;
    auto out_itr = std::begin(kmer_counts);
    auto kmer = 0;
    for (auto i = 0; i < kmer_size - 1 && seq_itr != seq_end; ++i)
      shift_kmer(kmer, *seq_itr++);
    while (seq_itr != seq_end) {
      shift_kmer(kmer, *seq_itr++);
      ++(*out_itr++)[kmer];
    }
    next_kmer_read += kmer_step;
  }

  auto
  operator+=(const kmer_counter &rhs) -> const kmer_counter &;

  [[nodiscard]] auto
  string(const std::uint64_t n_reads = 1) const -> std::string;

  [[nodiscard]] static auto
  decode_kmer(auto word, const auto n_bases);
};

template <> struct std::formatter<kmer_counter> : std::formatter<std::string> {
  auto
  format(const kmer_counter &kc, auto &ctx) const {
    return std::formatter<std::string>::format(kc.string(), ctx);
  }
};

#endif  // SRC_KMER_COUNTER_HPP_
