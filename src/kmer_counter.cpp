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

#include "kmer_counter.hpp"
#include "falco_utils.hpp"

#include <array>
#include <cstdint>
#include <vector>

auto
kmer_counter::operator+=(const kmer_counter &rhs) -> const kmer_counter & {
  two_dim_add(kmer_counts, rhs.kmer_counts);
  return *this;
}

[[nodiscard]] auto
kmer_counter::string(const std::uint64_t n_reads) const -> std::string {
  constexpr auto as_frac = [](const auto a, const auto b) {
    return static_cast<double>(a) / static_cast<double>(b);
  };
  auto r = std::format(">>Kmer Content\t{}\n", "pass");
  auto total_kmer_count = 0ul;
  for (const auto &c : kmer_counts)
    for (auto i = 0; i < n_kmers; ++i)
      total_kmer_count += c[i];
  r += std::format("total kmers: {}\n", total_kmer_count);
  return r + end_module_tag;
}

[[nodiscard]] auto
kmer_counter::decode_kmer(auto word, const auto n_bases) {
  static constexpr auto mask = 3u;
  static constexpr auto bits_per_base = 2u;
  static constexpr auto bases = "ACTG";
  std::string r;
  for (auto i = 0u; i < n_bases; ++i, word >>= bits_per_base)
    r += bases[word & mask];
  return r;
}
