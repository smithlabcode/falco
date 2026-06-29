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

#include "falco_utils.hpp"

#include "nlohmann/json.hpp"

#include <array>
#include <cmath>
#include <cstdint>
#include <span>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

[[nodiscard]] auto
file_info::string() const -> std::string {
  static constexpr auto n_indent = 4;
  nlohmann::json data = *this;
  return data.dump(n_indent);
}

[[nodiscard]] auto
size_to_units(const std::int64_t s) -> std::string {
  const auto as_frac_3 = [](const auto a, const auto b) {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
    return std::floor(100 * as_frac(a, b)) / 100;
  };
  if (s >= gigabytes)
    return std::format("{}GiB", as_frac_3(s, gigabytes));
  if (s >= megabytes)
    return std::format("{}MiB", as_frac_3(s, megabytes));
  if (s >= kilobytes)
    return std::format("{}KiB", as_frac_3(s, kilobytes));
  return std::format("{}", s);
};

[[nodiscard]] static auto
get_linear_interval(const std::uint64_t n_bases,
                    const std::uint64_t n_groups_target)
  -> std::tuple<std::uint64_t, std::uint64_t> {
  // returns the interval size and the corresponding number of groups
  static constexpr auto max_n_bases = 10'000'000LU;
  static constexpr auto multiplier = 10LU;
  static constexpr auto base_values = std::array{2LU, 5LU, 10LU};
  if (n_bases == 0 || n_bases > max_n_bases)
    throw std::runtime_error(
      std::format("bad n_bases in get_linear_interval: {}", n_bases));
  auto scaling_factor = 1;
  while (true) {
    for (const auto b : base_values) {
      const auto interval = b * scaling_factor;
      const auto n_groups = (n_bases + interval - 1) / interval;
      if (n_groups <= n_groups_target)
        return {interval, n_groups};
    }
    scaling_factor *= multiplier;
  }
  return {1LU, n_bases};
}

[[nodiscard]] auto
make_base_groups(const std::uint64_t n_bases, const std::uint64_t n_initial,
                 const std::uint64_t n_groups_target)
  -> std::vector<base_group_t> {
  static constexpr auto make_one_group = [](const auto a, const auto b) {
    return base_group_t{a, b};
  };
  static constexpr auto make_groups = [](const auto n, const auto offset,
                                         const auto scale) {
    const auto f = [scale, offset](const auto x) { return offset + x * scale; };
    return std::views::transform(std::views::iota(0LU, n + 1), f) |
           std::views::adjacent_transform<2>(make_one_group);
  };
  if (n_bases <= n_initial + n_groups_target)
    return make_groups(n_bases, 0, 1) | std::ranges::to<std::vector>();
  const auto n_non_initial = n_bases - n_initial;
  const auto [group_size, n_groups] =
    get_linear_interval(n_non_initial, n_groups_target);
  auto groups = make_groups(n_initial, 0, 1) | std::ranges::to<std::vector>();
  const auto groups1 = make_groups(n_groups, n_initial, group_size);
  groups.reserve(std::size(groups) + n_groups);
  std::ranges::copy(groups1, std::back_inserter(groups));
  // ADS: the line below is to try and make the same intervals as FastQC and
  // maybe isn't needed
  groups.back().second = n_bases;
  return groups;
}

[[nodiscard]] auto
get_default_base_groups(const std::uint64_t n_bases, const bool use_target)
  -> const std::vector<std::pair<std::uint64_t, std::uint64_t>> & {
  static constexpr auto default_n_initial = 9UL;
  static constexpr auto default_n_groups_target = 75UL - default_n_initial;
  static const auto base_group =
    use_target
      ? make_base_groups(n_bases, default_n_initial, default_n_groups_target)
      : make_base_groups(n_bases, 0UL, n_bases);
  return base_group;
}

[[nodiscard]] auto
make_group_tag(const base_group_t g) -> std::string {
  return (g.first + 1 == g.second)
           ? std::format("{}", g.first + 1)
           // ADS: make a closed interval
           : std::format("{}-{}", g.first + 1, g.second);
}

[[nodiscard]] auto
make_group_tag_quoted(const base_group_t g) -> std::string {
  return (g.first + 1 == g.second)
           ? std::format(R"({})", g.first + 1)
           // ADS: make a closed interval
           : std::format(R"({}-{})", g.first + 1, g.second);
}

[[nodiscard]] auto
get_theoretical_distribution(const falco::gc_content_array &gc,
                             const std::uint64_t total_count)
  -> std::vector<double> {
  static constexpr auto mode_width = 0.90;  // ADS: where does this come from?
  // calculate deviation of a hist from a "normal" with same mode and sd
  const auto n_bins = std::size(gc);
  const auto gc_beg = std::cbegin(gc);
  if (total_count <= 1)  // we will be dividing by (total_count - 1)
    return std::vector<double>(std::size(gc));

  // get mode
  const auto mode_itr = std::ranges::max_element(gc);
  const std::uint64_t mode_pos = std::distance(gc_beg, mode_itr);
  const auto mode_val = static_cast<double>(*mode_itr);

  // ADS: in case mode is not sharp average nearby values (not clear on why)
  const auto gt_cut = [&](const double x) { return x < mode_width * mode_val; };
  const auto right_itr = std::find_if(mode_itr, std::cend(gc), gt_cut);
  const auto left_itr =
    std::find_if(std::reverse_iterator(mode_itr), std::crend(gc), gt_cut);
  const auto n = std::distance(std::reverse_iterator(right_itr), left_itr);
  const auto mode = static_cast<double>(mode_pos) + as_frac(n - 1, 2.0);

  // theoretical distribution
  const auto cntr_sq = [m = mode](const auto v) { return (v - m) * (v - m); };
  const auto sd_term = [&](const auto &x) {
    return cntr_sq(std::get<0>(x)) * std::get<1>(x);
  };
  const auto id_gc = std::views::enumerate(gc) | std::views::transform(sd_term);
  const auto sd = std::sqrt(as_frac(
    std::reduce(std::cbegin(id_gc), std::cend(id_gc)), total_count - 1));
  const auto to_normal = [&](const auto val) {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
    return std::exp(-cntr_sq(val) / (2.0 * sd * sd));
  };
  auto normed = std::views::iota(0u, n_bins) |
                std::views::transform(to_normal) |
                std::ranges::to<std::vector>();
  const auto denom = std::reduce(std::cbegin(normed), std::cend(normed));
  std::ranges::transform(normed, std::begin(normed),
                         [&](const auto x) { return x * total_count / denom; });
  return normed;
}

[[nodiscard]] auto
sum_deviation_from_normal(const falco::gc_content_array &gc) -> double {
  const auto gc_beg = std::cbegin(gc);
  const auto gc_end = std::cend(gc);
  const auto total_count = std::reduce(gc_beg, gc_end);
  if (total_count <= 1)  // we will be dividing by (total_count - 1)
    return 0.0;

  const auto theor = get_theoretical_distribution(gc, total_count);
  const auto diff = [](const auto a, const auto b) { return std::fabs(a - b); };
  const auto r = std::transform_reduce(gc_beg, gc_end, std::cbegin(theor), 0.0,
                                       std::plus{}, diff);
  return pct(as_frac(r, total_count));
}

[[nodiscard]] auto
smooth_gc_content(const falco::gc_content_array &data,
                  const std::int64_t window_size) -> std::vector<double> {
  const auto get_mean = [&](std::ranges::viewable_range auto &&r) {
    return as_frac(std::reduce(std::cbegin(r), std::cend(r)), std::size(r));
  };
  assert(window_size < std::ssize(data));
  std::vector<double> smoothed;
  for (auto w = 1; w < (window_size + 1) / 2; ++w)
    smoothed.push_back(get_mean(
      std::ranges::subrange(std::cbegin(data), std::cbegin(data) + w)));
  for (const auto &window : data | std::views::slide(window_size))
    smoothed.push_back(get_mean(window));
  for (auto w = (window_size + 1) / 2; w > 1; --w)
    smoothed.push_back(get_mean(
      std::ranges::subrange(std::cend(data) - w + 1, std::cend(data))));
  return smoothed;
}
