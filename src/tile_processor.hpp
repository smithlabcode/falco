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

#include "bam_file.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"

#include <array>
#include <charconv>
#include <cstdint>
#include <iterator>
#include <limits>
#include <ranges>
#include <string>
#include <system_error>
#include <thread>  // IWYU pragma: keep
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

// ADS: notes
//
// - The number of tiles should never exceed 2500
// - Among the first 10k reads, all should contribute to tiles?
// - The quality score offset should be subtracted before any summary analysis
//   (e.g., computing the mean), but with enough data it doesn't seem to matter

struct tile_processor {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{5.0, "pass"},
    std::pair{10.0, "warn"},
    std::pair{std::numeric_limits<double>::max(), "fail"},
  };

  using qual_vec = std::vector<std::pair<std::uint64_t, std::uint64_t>>;
  static constexpr auto read_skip = 10 - 1;

  static std::uint32_t preceding_colons;

  std::int32_t read_idx{};
  std::uint32_t max_read_len{};
  std::uint32_t tile_id{};
  qual_vec::iterator qual{};
  std::unordered_map<std::uint32_t, qual_vec> quals;

  auto
  resize(const std::uint32_t updated_length) {
    max_read_len = updated_length;
    for (auto &q : quals | std::views::values)
      q.resize(max_read_len);
    const auto tile_id_itr = quals.find(tile_id);
    if (tile_id_itr != std::cend(quals))
      qual = std::begin(quals[tile_id]);  // resize invalidates iterators
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
        quals.emplace(tile_id, qual_vec(max_read_len, {0, 0}));
      qual = std::begin(quals[tile_id]);
    }
  }

  static auto
  set_preceding_colons(const std::string &fastq_filename) -> std::uint32_t;

  [[nodiscard]] auto
  string(const std::uint32_t len) const -> std::string;

  auto
  operator+=(const tile_processor &rhs) -> const tile_processor &;

  auto
  operator()(const auto &rec) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      const auto curr_len = static_cast<std::uint32_t>(get_seq_size(rec));
      if (curr_len > max_read_len)
        resize(curr_len);
      update_tile_id(get_name(rec), get_name_end(rec));
      if constexpr (std::is_same_v<std::decay_t<decltype(rec)>, bamrec>) {
        if (rec.is_rev)
          count_quals_itr_rev(get_qual(rec), get_qual_end(rec), qual);
        else
          count_quals_itr(get_qual(rec), get_qual_end(rec), qual);
      }
      else
        count_quals_itr(get_qual(rec), get_qual_end(rec), qual);
    }
  }
};

#endif  // SRC_TILE_PROCESSOR_HPP_
