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

#include "duplication_results.hpp"

#include "falco_utils.hpp"
#include "falco_word.hpp"

#include <algorithm>
#include <cstdint>
#include <format>
#include <functional>
#include <ranges>
#include <string>
#include <unordered_map>
#include <vector>

auto
duplication_results::operator+=(const duplication_results &rhs)
  -> const duplication_results & {
  for (const auto &[k, v] : rhs.dups)
    dups[k] += v;
  add_unique_seqs = add_unique_seqs && rhs.add_unique_seqs;
  n_unique = std::size(dups);
  limit_count += rhs.limit_count;
  return *this;
}

[[nodiscard]] auto
duplication_results::format_overrepresented(const std::uint64_t n_reads) const
  -> std::string {
  static constexpr auto start_tag = ">>Overrepresented sequences\t{}\n";
  static constexpr auto header = "#Sequence\t"
                                 "Count\t"
                                 "Percentage\t"
                                 "Possible Source\n";
  auto r = std::format(start_tag, "warn");
  r += header;
  const auto cutoff =
    static_cast<double>(n_reads) * duplication_results::overrepresented_cutoff;
  std::vector<std::pair<std::uint64_t, falco_word>> overrep_seqs;
  for (const auto &[seq, seq_count] : dups)
    if (seq_count >= cutoff)
      overrep_seqs.emplace_back(seq_count, seq);
  std::ranges::sort(overrep_seqs, std::greater{});
  for (const auto &[seq_count, seq] : overrep_seqs)
    r += std::format("{}\t{}\t{}\n", seq, seq_count,
                     pct(as_frac(seq_count, n_reads)));
  return r + end_module_tag;
}

[[nodiscard]] auto
duplication_results::format_duplication_levels() const -> std::string {
  static constexpr auto start_tag = ">>Sequence Duplication Levels\t{}\n"
                                    "#Total Deduplicated Percentage\t{}\n";
  static constexpr auto header = "#Duplication Level\t"
                                 "Percentage of deduplicated\t"
                                 "Percentage of total\n";
  // clang-format off
  static constexpr auto bin_breaks = std::array{
    1,
    2,
    3,
    4,
    5,
    6,
    7,
    8,
    9,
    10,
    50,
    100,
    500,
    1'000,
    5'000,
    10'000,
  };
  static constexpr auto bin_labels = std::array{
    "1",
    "2",
    "3",
    "4",
    "5",
    "6",
    "7",
    "8",
    "9",
    ">10",
    ">50",
    ">100",
    ">500",
    ">1k",
    ">5k",
    ">10k",
  };
  // clang-format on
  auto r = std::format(start_tag, "pass", 0.0);
  r += header;
  const auto max_count = std::ranges::max(std::views::values(dups));
  std::vector<std::uint64_t> hist(max_count + 1);
  for (const auto seq_count : std::views::values(dups))
    ++hist[seq_count];

  const auto seq_total = std::reduce(std::cbegin(hist), std::cend(hist));
  [[maybe_unused]] auto seq_dedup = 0ul;
  for (const auto [times_duped, seq_count] : std::views::enumerate(hist))
    seq_dedup += times_duped * seq_count;

  std::vector<std::uint64_t> binned(std::size(bin_breaks), 0);
  auto bin_id = 0ul;
  for (const auto [i, h] : std::views::enumerate(hist)) {
    if (bin_id + 1 < std::size(bin_breaks) && i > bin_breaks[bin_id])
      ++bin_id;
    binned[bin_id] += h;
  }

  for (const auto [lbl, bin] : std::views::zip(bin_labels, binned))
    r += std::format("{}\t{:.6g}\n", lbl, as_frac(bin, seq_total));

  return r;
}
