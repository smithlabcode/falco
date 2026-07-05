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

#ifndef SRC_ADAPTER_MATCHER_HPP_
#define SRC_ADAPTER_MATCHER_HPP_

#include "base_groups.hpp"
#include "falco_utils.hpp"

#include <cstdint>
#include <format>
#include <string>
#include <vector>

struct file_grades;
class run_mode;

struct adapter_matcher {
  std::uint32_t n_adapters{};
  std::uint32_t adapter_size{};
  std::vector<std::uint64_t> encoded_adapters;
  std::vector<std::vector<std::uint64_t>> adap_counts;

  auto
  apply_groups(const run_mode &mode) -> void;

  adapter_matcher();

  std::vector<std::uint64_t> encoded_read;
  std::vector<std::uint64_t>::const_iterator enc_beg;

  auto
  match_adapters(const auto seq, const auto len) {
    static const auto adap_mask = (1ul << adapter_size * nibble_size) - 1ul;
    if (len < adapter_size) [[unlikely]]
      return;
    if (std::size(encoded_read) < len) [[unlikely]] {
      encoded_read.resize(len);
      enc_beg = std::cbegin(encoded_read);
    }

    std::uint64_t enc{};
    auto i = 0u;
    while (i + 1 < adapter_size)
      enc = (enc << nibble_size) + encode_nibble(seq[i++]);
    auto idx = 0;
    while (i < len) {
      enc = (enc << nibble_size) + encode_nibble(seq[i++]);
      encoded_read[idx++] = (enc & adap_mask);
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
  operator+=(const adapter_matcher &rhs) -> const adapter_matcher &;

  [[nodiscard]] auto
  get_grade(const std::uint64_t n_reads) const -> std::string;

  [[nodiscard]] auto
  report(const std::uint64_t n_reads, const std::uint64_t max_read_len,
         const std::vector<base_group_t> &groups,
         const file_grades &grades) const -> std::string;

  [[nodiscard]] auto
  html(const std::uint64_t n_reads, const std::uint64_t max_read_len,
       const std::vector<base_group_t> &groups,
       const file_grades &grades) const -> std::string;
};

#endif  // SRC_ADAPTER_MATCHER_HPP_
