// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bamrec.hpp"

#include <algorithm>
#include <cstring>
#include <format>
#include <memory>
#include <span>
#include <string>
#include <utility>

[[nodiscard]] auto
bamrec::to_string() const -> std::string {
  const auto buffer_s = std::span(std::cbegin(buffer), std::cend(buffer));
  const auto name_itr = std::cbegin(buffer_s);
  const auto seq_itr = name_itr + name_len;
  const auto qual_itr = seq_itr + seq_len;
  std::string qual_fixed(seq_len, '\0');
  std::transform(qual_itr, qual_itr + seq_len, std::begin(qual_fixed),
                 [](const auto c) { return c + quality_score_offset; });
  return std::format("@{}\n{}\n+\n{}",                //
                     std::string(name_itr, seq_itr),  //
                     std::string(seq_itr, qual_itr),  //
                     qual_fixed);
}

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
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
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
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    *d_first++ = seq_nt16_str[bam_seqi(first, j)];
  return d_first;
}

[[nodiscard]] auto
bamrec::get_next(bamrec::pos_t &itr, const bamrec::pos_t end,
                 bamrec &rec) -> bool {
  if (std::distance(itr, end) < bam_core_t::sz)
    return false;
  bam_core_t core{};
  std::memcpy(std::addressof(core), std::to_address(itr), sizeof(bam_core_t));
  if (std::distance(itr, end) < core.real_block_size())
    return false;
  rec.name_len = core.l_read_name - 1;  // we don't need the '\0'
  rec.seq_len = core.l_seq;
  const auto rec_size = rec.name_len + 2 * rec.seq_len;
  if (std::size(rec.buffer) < rec_size)
    rec.buffer.resize(rec_size);
  auto out_itr = std::begin(rec.buffer);
  std::copy_n(itr + bam_core_t::read_name_offset, rec.name_len, out_itr);
  out_itr += rec.name_len;  // increment data cursor to sequence
  auto seq_in = itr + core.seq_offset();
  if (core.bam_is_rev())
    assign_sequence_revcomp(seq_in, core.l_seq, out_itr);
  else
    assign_sequence(seq_in, core.l_seq, out_itr);
  out_itr += rec.seq_len;  // increment data cursor to qual
  const auto has_qual = (itr[core.qual_offset()] != qual_missing_code);
  if (has_qual) {
    const auto qual_itr = itr + core.qual_offset();
    const auto qual_end = qual_itr + rec.seq_len;
    if (core.bam_is_rev())
      std::reverse_copy(qual_itr, qual_end, out_itr);
    else
      std::copy(qual_itr, qual_end, out_itr);
  }
  else
    *out_itr = static_cast<char>(qual_missing_code);
  itr += core.real_block_size();
  return true;
}

[[nodiscard]] auto
bamrec::find_end_pos(bamrec::pos_t itr,
                     const bamrec::pos_t end) -> bamrec::pos_t {
  static constexpr std::int64_t record_size_size = sizeof(std::uint32_t);
  std::uint32_t record_size{};
  while (itr != end) {
    if (std::distance(itr, end) < record_size_size)
      return itr;
    std::memcpy(std::addressof(record_size), std::to_address(itr),
                record_size_size);
    record_size += record_size_size;
    if (std::distance(itr, end) < record_size)
      return itr;
    itr += record_size;  // only increment on consume of full record
  }
  return itr;
}
