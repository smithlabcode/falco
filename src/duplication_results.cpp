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
  // overrepresented sequences
  auto r = std::format(">>Overrepresented sequences\t{}\n", "warn");
  // Keep only sequences that pass the input cutoff
  auto overrep_seqs = std::vector<std::pair<std::uint64_t, falco_word>>{};
  const auto overrep = static_cast<std::uint64_t>(
    n_reads * duplication_results::overrepresented_cutoff);
  for (const auto &[seq, seq_count] : dups)
    if (seq_count > overrep)
      overrep_seqs.emplace_back(seq_count, seq);
  std::ranges::sort(overrep_seqs, std::greater{});
  for (const auto &[seq_count, seq] : overrep_seqs)
    r += std::format("{}\t{}\t{}\n", seq, seq_count,
                     pct(as_frac(seq_count, n_reads)));
  return r + end_module_tag;
}

auto
duplication_results::summarize() {
  // Remove all but the top max_unique sequences
  const auto dups_sz = std::size(dups);
  if (dups_sz <= max_unique)
    return;
  std::vector<std::uint32_t> counter;
  counter.resize(dups_sz);
  for (const auto &d : dups)
    counter.push_back(d.second);
  std::ranges::sort(counter, std::ranges::greater{});
  const auto cutoff = counter[max_unique];
  auto dups_itr = std::cbegin(dups);
  while (dups_itr != std::cend(dups))
    if (dups_itr->second <= cutoff)
      dups_itr = dups.erase(dups_itr);
    else
      ++dups_itr;
  n_unique = std::size(dups);
}

[[nodiscard]] auto
duplication_results::string() const -> std::string {
  return std::format(">>Sequence Duplication Levels\t{}\n", "fail") +
         std::format("n_unique: {}\n"
                     "limit_count: {}\n",
                     n_unique, limit_count) +
         end_module_tag;
}
