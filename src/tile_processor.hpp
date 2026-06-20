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

#include "boost/boost_unordered.hpp"

#include <array>
#include <charconv>
#include <cstdint>
#include <format>
#include <iterator>
#include <limits>
#include <map>
#include <ranges>
#include <string>
#include <system_error>
#include <thread>  // IWYU pragma: keep
#include <type_traits>
#include <utility>
#include <vector>

namespace falco {
enum class encoding : std::uint8_t;
}

// Notes
//
// - The number of tiles should never exceed 2500?
// - Among the first 10k reads, all should contribute to tiles?

struct tile_processor {
  static constexpr auto grade_cutoffs = std::array{
    std::pair{5.0, "pass"},
    std::pair{10.0, "warn"},
    std::pair{std::numeric_limits<double>::max(), "fail"},
  };

  // ADS: needs to count roughly ~1M reads each contributing up to 128
  using qual_vec = std::vector<std::pair<std::uint64_t, std::uint64_t>>;
  static constexpr auto read_skip = 10 - 1;

  std::uint32_t tile_id_position{};
  std::int32_t read_idx{};
  std::uint32_t max_read_len{};
  std::uint32_t tile_id{};
  qual_vec::iterator qual{};
  boost::unordered_flat_map<std::uint32_t, qual_vec> quals;
  std::map<std::uint32_t, std::vector<double>> centered;

  auto
  trim() -> void;

  auto
  adjust_fastq_qual_encoding(const falco::encoding enc) -> void;

  /// finalize does 5 things:
  /// (1) trims tile data that's too long for the given tile
  /// (2) adjust the quality scores based on the encoding
  /// (3) makes groups for each tile's data vectors
  /// (4) calculates the (ordered) map of centered values
  /// (5) clears the 'quals' to force a crash if it's used afterwards
  auto
  finalize(const run_mode &mode, const file_info &info) -> void;

  auto
  init(const file_info &info) {
    tile_id_position = info.tile_id_position;
  }

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
    while (colon_count < tile_id_position && tile_itr != name_end)
      colon_count += (*tile_itr++ == ':');
    std::uint32_t curr_tile_id{};
    const auto [_, ec] = std::from_chars(tile_itr, name_end, curr_tile_id);
    if (ec != std::errc{})
      throw std::system_error(std::make_error_code(ec),
                              "failed to parse tile id");
    if (curr_tile_id != tile_id) {
      tile_id = curr_tile_id;
      auto tile_id_itr = quals.find(tile_id);
      if (tile_id_itr == std::cend(quals))
        tile_id_itr =
          quals.emplace(tile_id, qual_vec(max_read_len, {0, 0})).first;
      qual = std::begin(tile_id_itr->second);
    }
  }

  [[nodiscard]] auto
  get_grade() const -> std::string;

  [[nodiscard]] auto
  get_report(const std::vector<base_group_t> &groups,
             const std::string &grade) const -> std::string;

  [[nodiscard]] auto
  get_html(const std::vector<base_group_t> &groups) const -> std::string;

  auto
  operator+=(const tile_processor &rhs) -> const tile_processor &;

  auto
  operator()(const auto &rec) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      const auto curr_len = static_cast<std::uint32_t>(get_seq_size(rec));
      update_tile_id(get_name(rec), get_name_end(rec));
      if (curr_len > max_read_len)
        resize(curr_len);
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

auto
get_tile_info(const std::string &fastq_filename) -> std::uint32_t;

#endif  // SRC_TILE_PROCESSOR_HPP_
