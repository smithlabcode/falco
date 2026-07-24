// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BAMREC_HPP_
#define SRC_BAMREC_HPP_

#include <htslib/bgzf.h>   // for BGZF
#include <htslib/hfile.h>  // for htell

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <variant>
#include <vector>

#ifdef bam_is_rev
#undef bam_is_rev
#endif

#ifdef bam_seqi
#undef bam_seqi
#endif

#ifdef BAM_FREVERSE
#undef BAM_FREVERSE
#endif

class bamrec {
public:
  using pos_t = char *;

private:
  struct bamrec_core_t {
    static constexpr auto sz = 36;
    static constexpr auto read_name_offset = 36;
    std::uint32_t block_size{};  // 0
    std::int32_t refID{};        // 4
    std::int32_t pos{};          // 8
    std::uint8_t l_read_name{};  // 12
    std::uint8_t mapq{};         // 13
    std::uint16_t bin{};         // 14
    std::uint16_t n_cigar_op{};  // 16
    std::uint16_t flag{};        // 18
    std::uint32_t l_seq{};       // 20
    std::int32_t next_refID{};   // 24
    std::int32_t next_pos{};     // 28
    std::int32_t tlen{};         // 32
    // 36

    [[nodiscard]] constexpr auto
    seq_offset() const -> std::uint32_t {
      static constexpr auto cigar_op_size = sizeof(std::uint32_t);
      return read_name_offset + l_read_name + cigar_op_size * n_cigar_op;
    }

    [[nodiscard]] constexpr auto
    n_seq_bytes() const -> std::uint32_t {
      return (l_seq + 1) / 2;
    }

    [[nodiscard]] constexpr auto
    qual_offset() const -> std::uint32_t {
      return seq_offset() + n_seq_bytes();
    }

    [[nodiscard]] constexpr auto
    real_block_size() const -> std::uint32_t {
      return (sizeof block_size) + block_size;
    }

    [[nodiscard]] constexpr auto
    bam_is_rev() -> bool {
      static constexpr auto BAM_FREVERSE = 16;
      return (flag & BAM_FREVERSE) != 0;
    }
  };

  static constexpr auto quality_score_offset = 33;
  static constexpr auto qual_missing_code = 0xff;  // from sam.c
  static constexpr auto l_read_name_offset = 8;
  static constexpr auto n_cigar_op_offset = 12;
  static constexpr auto flag_offset = 14;
  static constexpr auto l_seq_offset = 16;

  bamrec_core_t core;
  std::vector<char> buffer;
  std::uint32_t name_len{};
  std::uint32_t seq_len{};

public:
  // clang-format off
  friend constexpr auto get_name(const bamrec &);
  friend constexpr auto get_name_end(const bamrec &);
  friend constexpr auto get_seq(const bamrec &);
  friend constexpr auto get_seq_end(const bamrec &);
  friend constexpr auto get_seq_size(const bamrec &);
  friend constexpr auto get_qual(const bamrec &);
  friend constexpr auto get_qual_end(const bamrec &);
  friend constexpr auto get_qual_size(const bamrec &);
  // clang-format on

  // clang-format off
  bamrec() = default;
  ~bamrec() = default;
  // delete copy and assignment
  bamrec(const bamrec &) = delete;
  auto operator=(const bamrec &) -> bamrec & = delete;
  auto operator=(bamrec &&) noexcept -> bamrec & = delete;
  // default move for emplace
  bamrec(bamrec &&) noexcept = delete; // default;
  // clang-format on

  [[nodiscard]] auto
  to_string() const {
    auto qual_itr = std::data(buffer) + name_len + seq_len;
    std::string qual_fixed(seq_len, '\0');
    std::transform(qual_itr, qual_itr + seq_len, std::begin(qual_fixed),
                   [](const auto c) { return c + quality_score_offset; });
    return std::format("@{}\n{}\n+\n{}", "name", "seq", qual_fixed);
  }

  operator bool() const { return name_len != 0; }

  [[nodiscard]] static auto
  get_next(auto &itr, const auto end, bamrec &rec) -> bool;

  [[nodiscard]] static auto
  find_end_pos(pos_t itr, const pos_t end) -> bamrec::pos_t;
};

[[nodiscard]] inline constexpr auto
get_name(const bamrec &rec) {
  return std::data(rec.buffer);
}

[[nodiscard]] inline constexpr auto
get_name_end(const bamrec &rec) {
  return get_name(rec) + rec.name_len;
}

[[nodiscard]] inline constexpr auto
get_seq(const bamrec &rec) {
  return get_name_end(rec);
}

[[nodiscard]] inline constexpr auto
get_seq_end(const bamrec &rec) {
  return get_seq(rec) + rec.seq_len;
}

[[nodiscard]] inline constexpr auto
get_seq_size(const bamrec &rec) {
  return rec.seq_len;
}

[[nodiscard]] inline constexpr auto
get_qual(const bamrec &rec) {
  return get_seq_end(rec);
}

[[nodiscard]] inline constexpr auto
get_qual_end(const bamrec &rec) {
  return get_qual(rec) + rec.seq_len;
}

[[nodiscard]] inline constexpr auto
get_qual_size(const bamrec &rec) {
  return get_seq_size(rec);
}

using bam_chunks_t = std::vector<std::pair<bamrec::pos_t, bamrec::pos_t>>;

struct bam_task_t {
  bamrec::pos_t beg{};
  bamrec::pos_t end{};
};

#ifdef seq_nt16_str
#undef seq_nt16_str
#endif

#ifdef bam_seqi
#undef bam_seqi
#endif

template <class BidirIt, class OutputIt>
static inline constexpr OutputIt
assign_sequence_revcomp(BidirIt first, auto last, OutputIt d_first) {
  constexpr auto complem = [](const auto x) {
    return "TNGNNNCNNNNNNNNNNNNA"[x - 'A'];
  };
  constexpr auto seq_nt16_str = "=ACMGRSVTWYHKDBN";
  constexpr auto bam_seqi = [](const auto s, const auto i) -> int {
    return s[i >> 1] >> ((~i & 1) << 2) & 0xf;
  };
  for (auto j = last; j != 0; ++d_first)
    *d_first = complem(seq_nt16_str[bam_seqi(first, --j)]);
  return d_first;
}

template <class BidirIt, class OutputIt>
static inline constexpr OutputIt
assign_sequence(BidirIt first, auto last, OutputIt d_first) {
  constexpr auto seq_nt16_str = "=ACMGRSVTWYHKDBN";
  constexpr auto bam_seqi = [](const auto s, const auto i) -> int {
    return s[i >> 1] >> ((~i & 1) << 2) & 0xf;
  };
  for (auto j = 0U; j != last; ++j)
    *d_first++ = seq_nt16_str[bam_seqi(first, j)];
  return d_first;
}

[[nodiscard]] inline auto
bamrec::get_next(auto &itr, const auto end, bamrec &rec) -> bool {
  bamrec::bamrec_core_t &core = rec.core;
  if (std::distance(itr, end) < bamrec::bamrec_core_t::sz)
    return false;
  std::memcpy(&core, itr, sizeof core);
  if (std::distance(itr, end) < core.real_block_size())
    return false;
  rec.name_len = (core.l_read_name - 1);  // we don't need the '\0'
  rec.seq_len = core.l_seq;
  const auto rec_size = rec.name_len + 2 * rec.seq_len;
  if (std::size(rec.buffer) < rec_size)
    rec.buffer.resize(rec_size);
  auto data_itr = std::data(rec.buffer);
  std::memcpy(data_itr, itr + bamrec_core_t::read_name_offset, rec.name_len);
  data_itr += rec.name_len;  // increment data cursor to sequence
  auto seq_in = itr + core.seq_offset();
  if (core.bam_is_rev())
    assign_sequence_revcomp(seq_in, core.l_seq, data_itr);
  else
    assign_sequence(seq_in, core.l_seq, data_itr);
  data_itr += rec.seq_len;  // increment data cursor to qual
  const auto has_qual = (itr[core.qual_offset()] != qual_missing_code);
  if (has_qual) {
    const auto qual_itr = itr + core.qual_offset();
    if (core.bam_is_rev())
      std::reverse_copy(qual_itr, qual_itr + rec.seq_len, data_itr);
    else
      // std::memcpy(data_itr, itr + core.qual_offset(), rec.seq_len);
      std::copy(qual_itr, qual_itr + rec.seq_len, data_itr);
  }
  else
    data_itr[0] = static_cast<char>(qual_missing_code);
  itr += core.real_block_size();
  return true;
}

[[nodiscard]] inline auto
get_next(auto &itr, const auto end, bamrec &rec) -> bool {
  return bamrec::get_next(itr, end, rec);
}

[[nodiscard]] inline auto
bamrec::find_end_pos(pos_t itr, const pos_t end) -> bamrec::pos_t {
  static constexpr std::int64_t record_size_size = sizeof(std::uint32_t);
  std::uint32_t record_size{};
  while (itr != end) {
    if (std::distance(itr, end) < record_size_size)
      return itr;
    std::memcpy(&record_size, itr, record_size_size);
    record_size += record_size_size;
    if (std::distance(itr, end) < record_size)
      return itr;
    itr += record_size;  // only increment on consume of full record
  }
  return itr;
}

#endif  // SRC_BAMREC_HPP_
