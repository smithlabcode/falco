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

#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <functional>
#include <numeric>
#include <ranges>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

[[nodiscard]] auto
kmer_result::operator<=>(const kmer_result &rhs) const {
  return obs_exp < rhs.obs_exp   ? std::strong_ordering::less
         : obs_exp > rhs.obs_exp ? std::strong_ordering::greater
                                 : std::strong_ordering::equal;
}

[[nodiscard]] auto
kmer_result::string() const {
  // ADS: format a row for output of a kmer
  return std::format("{}\t{}\t{:.6g}\t{:.6g}\t{}",
                     kmer_counter::decode_kmer(kmer, kmer_counter::kmer_size),
                     count, pval, obs_exp, pos + 1);
}

auto
kmer_counter::operator+=(const kmer_counter &rhs) -> const kmer_counter & {
  two_dim_add(kmer_counts, rhs.kmer_counts);
  return *this;
}

// ADS: series representation for the lower incomplete gamma P(a,x)
[[nodiscard]] static inline auto
gamma_p_series(const double a, const double x) -> double {
  static constexpr auto eps = std::numeric_limits<double>::epsilon();
  static constexpr auto max_iter = 100;
  double sum = 1.0 / a;
  double term = sum;
  for (auto n = 1; n < max_iter; ++n) {
    term *= x / (a + n);
    sum += term;
    if (term < eps * sum)
      break;
  }
  return sum * std::exp(-x + a * std::log(x) - std::lgamma(a));
}

[[nodiscard]] static constexpr inline auto
safe_floor(const auto x, const auto floor_val) {
  return std::abs(x) < floor_val ? floor_val : x;
}

// ADS: continued fraction representation for the upper incomplete gamma Q(a,x)
[[nodiscard]] static inline auto
gamma_q_contfrac(const double a, const double x) {
  // ADS: fpmin is 2.2e-16 above smallest positive double
  static constexpr auto epsilon = std::numeric_limits<double>::epsilon();
  static constexpr auto fpmin = std::numeric_limits<double>::min() / epsilon;
  static constexpr auto max_iter = 100;
  double b = x + 1.0 - a;
  double c = 1.0 / fpmin;
  double d = 1.0 / b;
  double h = d;
  for (auto i = 1; i < max_iter; ++i) {
    const double an = -i * (i - a);
    b += 2;
    d = safe_floor(an * d + b, fpmin);
    c = safe_floor(b + an / c, fpmin);
    d = 1.0 / d;
    const double delta = d * c;
    h *= delta;
    if (std::abs(delta - 1.0) < epsilon)
      break;
  }
  return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
}

// ADS: regularized lower incomplete gamma P(a,x)
[[nodiscard]] static inline auto
gamma_p(const double a, const double x) {
  return x <= 0.0 || a <= 0 ? 0.0
                            : (x < a + 1.0 ? gamma_p_series(a, x)
                                           : 1.0 - gamma_q_contfrac(a, x));
}

[[nodiscard]] static inline auto
poisson_q(const std::uint64_t k, const double mu) {
  return mu > 0.0 ? gamma_p(static_cast<double>(k) + 1.0, mu) : 1.0;
}

[[nodiscard]] auto
kmer_counter::get_kmer_results() const -> std::vector<kmer_result> {
  const auto range_reduce = [](const auto &v) {
    return std::reduce(std::cbegin(v), std::cend(v));
  };

  auto counts_by_kmer = std::vector<std::uint64_t>(n_kmers);
  std::ranges::for_each(kmer_counts,
                        [&](const auto &k) { add(counts_by_kmer, k); });
  auto counts_by_pos = std::vector<std::uint64_t>(max_read_len);
  std::ranges::transform(kmer_counts, std::begin(counts_by_pos),
                         [&](const auto &a) { return range_reduce(a); });
  const auto total = static_cast<double>(range_reduce(counts_by_kmer));
  const auto to_probs = [&](const auto &v) {
    const auto p = [&](const auto x) { return as_frac(x, total); };
    return v | std::views::transform(p) | std::ranges::to<std::vector>();
  };

  const auto kmer_probs = to_probs(counts_by_kmer);
  const auto pos_probs = to_probs(counts_by_pos);

  const auto kr_maker = [kmer = 0ul](const auto c) mutable {
    return kmer_result{.kmer = kmer++, .count = c};
  };
  auto results = counts_by_kmer | std::views::transform(kr_maker) |
                 std::ranges::to<std::vector>();

  for (const auto [pos, pos_counts] : std::views::enumerate(kmer_counts)) {
    const auto count_for_pos = static_cast<double>(counts_by_pos[pos]);
    const auto pos_prob = pos_probs[pos];
    const auto per_kmer = std::views::zip(pos_counts, kmer_probs, results);
    for (auto [count, prob, res] : per_kmer) {
      const auto expected = total * prob * pos_prob;
      if (const auto obs_exp = as_frac(count, expected);
          obs_exp > res.obs_exp) {
        res.obs_exp = obs_exp;
        res.pos = pos;
      }
      const auto pval = poisson_q(count, prob * count_for_pos) * n_kmers;
      res.pval = pval < res.pval ? pval : res.pval;
    }
  }

  const auto to_erase = std::ranges::remove_if(results, [&](const auto &x) {
    return x.obs_exp < min_obs_exp_to_report || x.pval > max_pval_to_report;
  });
  results.erase(std::cbegin(to_erase), std::cend(to_erase));
  std::ranges::sort(results, std::greater{});
  return results;
}

[[nodiscard]] auto
kmer_counter::get_report(std::string &grade) const -> std::string {
  static constexpr auto start_tag = ">>Kmer Content\t{}\n";
  static constexpr auto header = "#Sequence\t"
                                 "Count\t"
                                 "PValue\t"
                                 "Obs/Exp Max\t"
                                 "Max Obs/Exp Position\n";
  const auto results = get_kmer_results();
  const auto pval = results.front().pval;
  const auto neg_log_p_val =
    pval > 0.0 ? -std::log10(pval) : std::numeric_limits<double>::max();
  grade = get_grade(grade_cutoffs, neg_log_p_val);
  auto r = std::format(start_tag, grade);
  r += header;
  for (const auto &res : results | std::views::take(n_kmers_to_report))
    // cppcheck-suppress useStlAlgorithm
    r += std::format("{}\n", res.string());
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
    // NOLINTNEXTLINE (cppcoreguidelines-pro-bounds-pointer-arithmetic)
    r += bases[word & mask];
  std::ranges::reverse(r);
  return r;
}
