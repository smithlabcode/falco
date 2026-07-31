// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_ADAPTER_MATCHER_HPP_
#define SRC_ADAPTER_MATCHER_HPP_

#include "base_groups.hpp"
#include "falco_utils.hpp"

#include <algorithm>
#include <cstdint>
#include <iterator>
#include <string>
#include <vector>
struct file_grades;

struct adapter_matcher {
  std::uint32_t n_adapters{};
  std::uint32_t adapter_size{};
  std::vector<std::uint64_t> encoded_adapters;
  std::vector<std::vector<std::uint64_t>> adap_counts;
  std::vector<std::uint64_t> encoded_read;

  auto
  apply_groups(const base_group_vec &groups) -> void;

  adapter_matcher();

  auto
  match_adapters(const auto seq, const auto len) {
    static const auto adap_mask = (1ul << adapter_size * nibble_size) - 1ul;
    if (len < adapter_size) [[unlikely]]
      return;
    if (std::size(encoded_read) < len) [[unlikely]]
      encoded_read.resize(len);

    const auto enc_beg = std::begin(encoded_read);
    std::uint64_t enc{};
    auto i = 0u;
    while (i + 1 < adapter_size)
      enc = (enc << nibble_size) + encode_nibble(seq[i++]);
    auto itr = enc_beg;
    while (i < len) {
      enc = (enc << nibble_size) + encode_nibble(seq[i++]);
      *itr++ = (enc & adap_mask);
    }
    const auto lim = enc_beg + len - adapter_size + 1;
    for (auto j = 0u; j < n_adapters; ++j)
      if (const auto p = std::find(enc_beg, lim, encoded_adapters[j]); p != lim)
        ++adap_counts[std::distance(enc_beg, p)][j];
  }

  auto
  resize(const std::uint32_t updated_length) {
    adap_counts.resize(updated_length,
                       std::vector<std::uint64_t>(n_adapters, 0));
  }

  auto
  release() {
    adap_counts.clear();
    adap_counts.shrink_to_fit();
  }

  auto
  operator+=(const adapter_matcher &rhs) -> const adapter_matcher &;

  auto
  add_and_consume(adapter_matcher &&rhs) -> void;

  [[nodiscard]] auto
  get_grade(const std::uint64_t n_reads) const -> std::string;

  [[nodiscard]] auto
  report(const std::uint64_t n_reads, const std::uint64_t max_read_len,
         const base_group_vec &groups,
         const file_grades &grades) const -> std::string;

  [[nodiscard]] auto
  html(const std::uint64_t n_reads, const std::uint64_t max_read_len,
       const base_group_vec &groups,
       const file_grades &grades) const -> std::string;
};

#endif  // SRC_ADAPTER_MATCHER_HPP_
