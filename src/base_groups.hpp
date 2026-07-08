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

#ifndef SRC_BASE_GROUPS_HPP_
#define SRC_BASE_GROUPS_HPP_

#include "falco_utils.hpp"

#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <iterator>
#include <ranges>
#include <string>
#include <utility>
#include <vector>

using base_group_t = std::pair<std::uint64_t, std::uint64_t>;
using base_group_vec = std::vector<base_group_t>;

[[nodiscard]] auto
make_base_groups(const std::uint64_t n_bases, const std::uint64_t n_initial,
                 const std::uint64_t n_groups_target) -> base_group_vec;

[[nodiscard]] auto
get_default_base_groups(const std::uint64_t n_bases,
                        const bool make_groups) -> base_group_vec;

[[nodiscard]] auto
make_group_tag(const base_group_t g) -> std::string;

[[nodiscard]] auto
make_group_tag_quoted(const base_group_t g) -> std::string;

[[nodiscard]] auto
apply_base_groups(const base_group_vec &groups, auto &rows) {
  assert(std::size(rows) <= groups.back().second);
  auto group_itr = std::cbegin(groups);
  auto current_row = 0U;
  for (const auto [idx, row] :
       std::views::enumerate(rows) | std::views::drop(1)) {
    if (static_cast<std::uint64_t>(idx) < group_itr->second)
      add(rows[current_row], row);
    else {
      std::swap(row, rows[++current_row]);
      ++group_itr;
    }
  }
  ++current_row;  // move past the last row used
  rows.resize(current_row);
}

[[nodiscard]] auto
apply_base_groups(const base_group_vec &groups, auto &rows, const auto &adder) {
  assert(std::size(rows) <= groups.back().second);
  auto group_itr = std::cbegin(groups);
  auto current_row = 0U;
  for (const auto [idx, row] :
       std::views::enumerate(rows) | std::views::drop(1)) {
    if (static_cast<std::uint64_t>(idx) < group_itr->second)
      adder(rows[current_row], row);
    else {
      std::swap(row, rows[++current_row]);
      ++group_itr;
    }
  }
  ++current_row;  // move past the last row used
  rows.resize(current_row);
}

#endif  // SRC_BASE_GROUPS_HPP_
