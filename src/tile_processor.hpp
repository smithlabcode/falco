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

#ifndef SRC_TILE_PROCESSOR_HPP_
#define SRC_TILE_PROCESSOR_HPP_

#include <charconv>
#include <cstdint>
#include <format>
#include <iterator>
#include <ranges>
#include <string>
#include <system_error>
#include <unordered_map>
#include <utility>
#include <vector>

struct tile_processor {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{0.05, "pass"},
    std::pair{0.10, "warn"},
    std::pair{1.00, "error"},
  };

  using qual_vec = std::vector<std::pair<std::uint64_t, std::uint64_t>>;

  static constexpr auto header = "#Tile\tBase\tMean\n";
  static constexpr auto tile_step = 10;
  static std::uint32_t preceding_colons;

  std::uint32_t next_tile_read{};
  std::uint32_t tile_id{};
  qual_vec::iterator qual{};
  std::unordered_map<std::uint32_t, qual_vec> quals;
  std::uint32_t max_read_len{};

  explicit tile_processor(const std::uint32_t max_read_len) :
    max_read_len{max_read_len} {}

  auto
  resize(const std::uint32_t updated_length) {
    max_read_len = updated_length;
    for (auto &q : quals | std::views::values)
      q.resize(max_read_len);
  }

  auto
  update_tile_id(const auto name_beg, const auto name_end) {
    auto tile_itr = name_beg;
    auto colon_count = 0u;
    while (colon_count < preceding_colons && tile_itr != name_end)
      colon_count += (*tile_itr++ == ':');
    std::uint32_t curr_tile_id{};
    const auto [_, ec] = std::from_chars(tile_itr, name_end, curr_tile_id);
    if (ec != std::errc{})
      throw std::system_error(std::make_error_code(ec),
                              "failed to parse tile id");
    if (curr_tile_id != tile_id) {
      tile_id = curr_tile_id;
      if (!quals.contains(tile_id))
        quals[tile_id] = qual_vec(max_read_len, {0, 0});
      qual = std::begin(quals[tile_id]);
    }
  }

  static auto
  set_preceding_colons(const std::string &fastq_filename) -> std::uint32_t;

  [[nodiscard]] auto
  string(const std::uint32_t len = 0) const -> std::string;

  auto
  operator+=(const tile_processor &rhs) -> const tile_processor &;

  auto
  operator()(auto q_itr, const auto q_end) {
    auto tab_itr = qual;
    while (q_itr != q_end) {
      tab_itr->first += *q_itr++;
      ++(*tab_itr++).second;
    }
  }
};

template <>
struct std::formatter<tile_processor> : std::formatter<std::string> {
  auto
  format(const tile_processor &tp, auto &ctx) const {
    return std::formatter<std::string>::format(tp.string(), ctx);
  }
};

#endif  // SRC_TILE_PROCESSOR_HPP_
