// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_TILE_PROCESSOR_HPP_
#define SRC_TILE_PROCESSOR_HPP_

#include "base_groups.hpp"
#include "falco_utils.hpp"
#include "quality_score.hpp"

#include "boost/boost_unordered.hpp"

#include <charconv>
#include <cstdint>
#include <iterator>
#include <map>
#include <ranges>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

struct file_info;

// Notes
//
// - The number of tiles should never exceed 2500?
// - Among the first 10k reads, all should contribute to tiles?

class tile_processor {
public:
  using tiles_centered_t = std::map<std::uint32_t, std::vector<double>>;
  // ADS: needs to count roughly ~1M reads each contributing up to 128
  using qual_vec = std::vector<std::pair<std::uint64_t, std::uint64_t>>;
  static constexpr auto read_skip = 10 - 1;

  std::uint32_t tile_id_position{};
  std::int32_t read_idx{};
  std::uint32_t max_read_len{};
  std::uint32_t tile_id{};
  qual_vec::iterator qual{};  // because we often don't change tile id
  boost::unordered_flat_map<std::uint32_t, qual_vec> quals;

  auto
  init(const file_info &info) -> void;

  /// finalize does 2 things:
  /// (1) trims tile data that's too long for the given tile
  /// (2) adjust the quality scores based on the encoding
  auto
  finalize(const file_info &info) -> void;

  [[nodiscard]] auto
  get_centered() -> tile_processor::tiles_centered_t;

  auto
  apply_groups(const base_group_vec &groups) -> void;

  auto
  operator()(const auto &rec) {
    if (read_idx-- == 0) [[unlikely]] {
      read_idx = read_skip;
      const auto curr_len = static_cast<std::uint32_t>(get_seq_size(rec));
      update_tile_id(get_name(rec), get_name_end(rec));
      if (curr_len > max_read_len)
        resize(curr_len);
      count_quals_itr(get_qual(rec), get_qual_end(rec), qual);
    }
  }

  auto
  operator+=(const tile_processor &rhs) -> const tile_processor &;

  auto
  add_and_consume(tile_processor &&rhs) -> void;

private:
  auto
  trim() -> void;

  auto
  adjust_fastq_qual_encoding(const file_info &info) -> void;

  auto
  resize(const std::uint32_t updated_length) -> void {
    max_read_len = updated_length;
    for (auto &q : quals | std::views::values)
      q.resize(max_read_len);
    const auto tile_id_itr = quals.find(tile_id);
    if (tile_id_itr != std::cend(quals))
      qual = std::begin(quals[tile_id]);  // resize invalidates iterators
  }

  auto
  release() -> void {
    quals.clear();
    boost::unordered_flat_map<std::uint32_t, qual_vec>().swap(quals);
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
};

auto
get_tile_info(const std::string &fastq_filename) -> std::uint32_t;

[[nodiscard]] auto
get_grade_tile(const tile_processor::tiles_centered_t &centered) -> std::string;

#endif  // SRC_TILE_PROCESSOR_HPP_
