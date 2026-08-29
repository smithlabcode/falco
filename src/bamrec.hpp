// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BAMREC_HPP_
#define SRC_BAMREC_HPP_

#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

#ifdef bam_is_rev
#undef bam_is_rev
#endif

#ifdef BAM_FREVERSE
#undef BAM_FREVERSE
#endif

struct bam_core_t {
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
};  // bam_core_t

class bamrec {
public:
  using pos_t = std::vector<char>::const_iterator;

private:
  static constexpr auto quality_score_offset = 33;
  static constexpr char qual_missing_code = -1;  // 0xff from sam.c

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
  get_next(pos_t &itr, const pos_t end, bamrec &rec) -> bool;

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

[[nodiscard]] inline auto
get_next(bamrec::pos_t &itr, const bamrec::pos_t end, bamrec &rec) -> bool {
  return bamrec::get_next(itr, end, rec);
}

#endif  // SRC_BAMREC_HPP_
