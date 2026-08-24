// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BAMREC_HPP_
#define SRC_BAMREC_HPP_

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iterator>
#include <span>
#include <string>
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
  using pos_t = std::vector<char>::const_iterator;

private:
  struct core_t {
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
  static constexpr std::uint8_t qual_missing_code = 0xff;  // from sam.c
  static constexpr auto l_read_name_offset = 8;
  static constexpr auto n_cigar_op_offset = 12;
  static constexpr auto flag_offset = 14;
  static constexpr auto l_seq_offset = 16;

  core_t core;
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
  bamrec(const bamrec &) = delete;
  auto operator=(const bamrec &) -> bamrec & = delete;
  auto operator=(bamrec &&) noexcept -> bamrec & = delete;
  bamrec(bamrec &&) noexcept = delete;
  // clang-format on

  [[nodiscard]] auto
  to_string() const -> std::string;

  operator bool() const { return name_len != 0; }

  [[nodiscard]] static auto
  get_next(auto &itr, const auto end, bamrec &rec) -> bool;

  [[nodiscard]] static auto
  find_end_pos(pos_t itr, const pos_t end) -> pos_t;
};

[[nodiscard]] inline constexpr auto
get_name(const bamrec &rec) {
  return std::cbegin(rec.buffer);
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
    constexpr auto low_nibble_on = 0xf;
    return s[i >> 1] >> ((~i & 1) << 2) & low_nibble_on;
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
    constexpr auto low_nibble_on = 0xf;
    return s[i >> 1] >> ((~i & 1) << 2) & low_nibble_on;
  };
  for (auto j = 0U; j != last; ++j)
    *d_first++ = seq_nt16_str[bam_seqi(first, j)];
  return d_first;
}

[[nodiscard]] inline auto
bamrec::get_next(auto &itr, const auto end, bamrec &rec) -> bool {
  if (std::distance(itr, end) < core_t::sz)
    return false;
  core_t &c = rec.core;
  std::memcpy(&c, std::data(std::span{itr, end}), sizeof c);
  if (std::distance(itr, end) < c.real_block_size())
    return false;
  rec.name_len = c.l_read_name - 1;  // we don't need the '\0'
  rec.seq_len = c.l_seq;
  const auto rec_size = rec.name_len + 2 * rec.seq_len;
  if (std::size(rec.buffer) < rec_size)
    rec.buffer.resize(rec_size);
  auto out_itr = std::begin(rec.buffer);
  std::copy_n(itr + core_t::read_name_offset, rec.name_len, out_itr);
  out_itr += rec.name_len;  // increment data cursor to sequence
  auto seq_in = itr + c.seq_offset();
  if (c.bam_is_rev())
    assign_sequence_revcomp(seq_in, c.l_seq, out_itr);
  else
    assign_sequence(seq_in, c.l_seq, out_itr);
  out_itr += rec.seq_len;  // increment data cursor to qual
  const auto has_qual = (itr[c.qual_offset()] != qual_missing_code);
  if (has_qual) {
    const auto qual_itr = itr + c.qual_offset();
    if (c.bam_is_rev())
      std::reverse_copy(qual_itr, qual_itr + rec.seq_len, out_itr);
    else
      std::copy(qual_itr, qual_itr + rec.seq_len, out_itr);
  }
  else
    *out_itr = static_cast<char>(qual_missing_code);
  itr += c.real_block_size();
  return true;
}

[[nodiscard]] inline auto
get_next(auto &itr, const auto end, bamrec &rec) -> bool {
  return bamrec::get_next(itr, end, rec);
}

[[nodiscard]] inline auto
bamrec::find_end_pos(pos_t itr, const pos_t end) -> pos_t {
  static constexpr std::int64_t record_size_size = sizeof(std::uint32_t);
  std::uint32_t record_size{};
  while (itr != end) {
    if (std::distance(itr, end) < record_size_size)
      return itr;
    std::memcpy(&record_size, std::data(std::span{itr, end}), record_size_size);
    record_size += record_size_size;
    if (std::distance(itr, end) < record_size)
      return itr;
    itr += record_size;  // only increment on consume of full record
  }
  return itr;
}

#endif  // SRC_BAMREC_HPP_
