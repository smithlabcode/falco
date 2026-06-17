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

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <tuple>
#include <vector>

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

[[nodiscard]] static auto
make_linear_base_groups(const std::uint64_t n_bases)
  -> std::vector<base_group_t> {
  static constexpr auto n_initial = 10LU;
  static constexpr auto n_groups_target = 65LU;
  static constexpr auto make_one_group = [](const auto a, const auto b) {
    return base_group_t{a, b};
  };
  static constexpr auto make_groups = [](const auto n, const auto offset,
                                         const auto scale) {
    return std::views::transform(
             std::views::iota(0LU, n + 1),
             [offset, scale](const auto x) { return offset + x * scale; }) |
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
  return groups;
}

[[nodiscard]] auto
make_base_groups(const std::int64_t n_bases)
  -> std::vector<std::pair<std::uint64_t, std::uint64_t>> {
  return make_linear_base_groups(n_bases);
}
