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
#include <compare>
#include <cstdint>
#include <format>
#include <string>
#include <vector>

auto
kmer_counter::operator+=(const kmer_counter &rhs) -> const kmer_counter & {
  two_dim_add(kmer_counts, rhs.kmer_counts);
  return *this;
}

struct kmer_result {
  std::uint64_t kmer{};
  std::uint64_t count{};
  double pval{};
  double obs_exp{};
  std::uint64_t pos{};
  [[nodiscard]] auto
  operator<=>(const kmer_result &rhs) const {
    return obs_exp < rhs.obs_exp   ? std::strong_ordering::less
           : obs_exp > rhs.obs_exp ? std::strong_ordering::greater
                                   : std::strong_ordering::equal;
  }
  [[nodiscard]] auto
  string() const {
    return std::format("{}\t{}\t{}\t{}\t{}",
                       kmer_counter::decode_kmer(kmer, kmer_counter::kmer_size),
                       count, pval, obs_exp, pos);
  }
};

[[nodiscard]] auto
kmer_counter::string([[maybe_unused]] const std::uint64_t n_reads) const
  -> std::string {
  static constexpr auto start_tag = ">>Kmer Content\t{}\n";
  static constexpr auto header = "#Sequence\t"
                                 "Count\t"
                                 "PValue\t"
                                 "Obs/Exp Max\t"
                                 "Max Obs/Exp Position\n";
  [[maybe_unused]] static constexpr auto max_kmers_to_plot = 10;
  static constexpr auto n_kmers_to_report = 20;
  static constexpr auto min_obs_exp_to_report = 5.0;

  auto counts_by_kmer = std::vector<std::uint64_t>(n_kmers);
  std::ranges::for_each(kmer_counts,
                        [&](const auto &k) { add(counts_by_kmer, k); });

  auto count_by_pos = std::vector<std::uint64_t>(max_read_len);
  std::ranges::transform(
    kmer_counts, std::begin(count_by_pos),
    [](const auto &x) { return std::reduce(std::cbegin(x), std::cend(x)); });

  std::vector<kmer_result> kmer_info(n_kmers);
  for (const auto [kmer, count] : std::views::enumerate(counts_by_kmer)) {
    kmer_info[kmer].kmer = kmer;
    kmer_info[kmer].count = count;
  }

  for (const auto [pos, pos_counts] : std::views::enumerate(kmer_counts))
    for (const auto [ki, kmer_count] : std::views::zip(kmer_info, pos_counts)) {
      // expected is expected for the current position
      const auto expected = as_frac(count_by_pos[pos], n_kmers);
      if (const auto oe = as_frac(kmer_count, expected); oe > ki.obs_exp) {
        ki.obs_exp = oe;
        ki.pos = pos;
      }
    }

  const auto to_erase = std::ranges::remove_if(kmer_info, [&](const auto &x) {
    return x.obs_exp < min_obs_exp_to_report;
  });
  kmer_info.erase(std::cbegin(to_erase), std::cend(to_erase));
  std::ranges::sort(kmer_info, std::greater{});

  const auto max_obs_exp = kmer_info.empty() ? 0.0 : kmer_info.front().obs_exp;
  auto r = std::format(start_tag, get_grade(grade_cutoffs, max_obs_exp));
  r += header;
  for (const auto &info : kmer_info | std::views::take(n_kmers_to_report))
    r += std::format("{}\n", info.string());
  return r + end_module_tag;
}

[[nodiscard]] auto
kmer_counter::decode_kmer(auto word, const auto n_bases) -> std::string {
  static constexpr auto mask = 3u;
  static constexpr auto bits_per_base = 2u;
  static constexpr auto bases = "ACTG";
  std::string r;
  using pos_t = std::decay_t<decltype(n_bases)>;
  for (pos_t i = 0; i < n_bases; ++i, word >>= bits_per_base)
    r += bases[word & mask];
  std::ranges::reverse(r);
  return r;
}
