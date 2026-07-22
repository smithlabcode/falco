// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BAMREC_HPP_
#define SRC_BAMREC_HPP_

#include <htslib/sam.h>

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <concepts>
#include <cstdint>
#include <cstdio>  // IWYU pragma: keep
#include <cstdlib>
#include <cstring>  // IWYU pragma: keep
#include <filesystem>
#include <format>  // IWYU pragma: keep
#include <iterator>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

static constexpr auto qual_missing_code = 0xff;  // from sam.c

struct bam_seq_itr {
  // ADS: this class needs to go
  const std::uint8_t *p{};
  std::uint64_t i{};
  explicit bam_seq_itr(const auto p) : p{p} {}
  [[nodiscard]] auto
  operator<=>(const bam_seq_itr &) const = default;
  [[nodiscard]] auto
  operator*() const {
    return seq_nt16_str[bam_seqi(p, i)];
  }
  [[nodiscard]] auto
  operator[](const auto j) const {
    return seq_nt16_str[bam_seqi(p, i + j)];
  }
  auto
  operator+(const auto rhs) {
    bam_seq_itr tmp(*this);
    tmp.i += rhs;
    return tmp;
  }
  auto
  operator++(int) {
    bam_seq_itr tmp(*this);
    ++i;
    return tmp;
  }
  auto
  operator++() {
    ++i;
    return *this;
  }
};

struct bamrec {
  using pos_t = std::vector<bam1_t>::const_iterator;
  std::int32_t l_qname{};
  std::int32_t l_qseq{};
  bool is_rev{};
  const char *n{};
  const std::uint8_t *r{};
  const std::uint8_t *q{};
  bamrec() = default;
  explicit bamrec(const bam1_t &b) :
    l_qname{b.core.l_qname}, l_qseq{b.core.l_qseq}, is_rev{bam_is_rev(&b)},
    n{bam_get_qname(&b)}, r{bam_get_seq(&b)}, q{bam_get_qual(&b)} {}
  [[nodiscard]] operator bool() const { return n != nullptr; }
};

// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
[[nodiscard]] constexpr auto
get_name(const bamrec &rec) {
  return rec.n;
}

[[nodiscard]] constexpr auto
get_name_end(const bamrec &rec) {
  return rec.n + rec.l_qname;
}

[[nodiscard]] constexpr auto
get_seq(const bamrec &rec) {
  return bam_seq_itr(rec.r);
}

[[nodiscard]] constexpr auto
get_seq_end(const bamrec &rec) {
  return bam_seq_itr(rec.r) + rec.l_qseq;
}

[[nodiscard]] constexpr auto
get_seq_size(const bamrec &rec) {
  return rec.l_qseq;
}

[[nodiscard]] constexpr auto
get_qual(const bamrec &rec) {
  return rec.q;
}

[[nodiscard]] constexpr auto
get_qual_end(const bamrec &rec) {
  return *rec.q == qual_missing_code ? rec.q : rec.q + rec.l_qseq;
}

[[nodiscard]] constexpr auto
get_qual_size(const bamrec &rec) {
  return rec.l_qseq;
}
// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

[[nodiscard]] inline auto
to_string(const bamrec &b) {
  // converts the bam record to FASTQ format for visualization
  static constexpr auto quality_score_offset = 33;
  static constexpr auto other = "\n+\n";
  auto name = std::format("@{}\n", std::string(b.n, b.l_qname));
  std::string read;
  auto seq_itr = get_seq(b);
  while (seq_itr != get_seq_end(b))
    read += *seq_itr++;
  std::string qual;
  auto qual_itr = get_qual(b);
  if (*qual_itr == qual_missing_code)
    qual = std::string(get_qual_size(b), 'B');
  else
    while (qual_itr != get_qual_end(b))
      qual += static_cast<char>(quality_score_offset + *qual_itr++);
  return name + read + other + qual + '\n';
}

struct bam_task_t {
  bamrec::pos_t beg{};
  bamrec::pos_t end{};
};

using bam_chunks_t = std::vector<std::pair<bamrec::pos_t, bamrec::pos_t>>;

[[nodiscard]] inline auto
get_next(bamrec::pos_t &cursor,
         [[maybe_unused]] const bamrec::pos_t end_itr) -> bamrec {
  auto tmp = cursor;
  ++cursor;
  return bamrec(*tmp);
}

#endif  // SRC_BAMREC_HPP_
