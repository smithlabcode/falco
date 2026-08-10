// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "samrec.hpp"

#include <algorithm>
#include <charconv>
#include <iterator>
#include <system_error>

#ifdef BAM_FREVERSE
#undef BAM_FREVERSE
#endif

#ifdef bam_is_rev
#undef bam_is_rev
#endif

[[nodiscard]] auto
samrec::get_next(samrec::pos_t &cursor, const samrec::pos_t end_itr,
                 samrec &rec) -> bool {
  static constexpr auto n_fields_to_skip = 7;
  constexpr auto complem = [](const auto x) {
    return "TNGNNNCNNNNNNNNNNNNA"[x - 'A'];
  };
  static constexpr auto bam_is_rev = [](const auto flag) -> bool {
    static constexpr auto BAM_FREVERSE = 16;
    return (flag & BAM_FREVERSE) != 0;
  };

  auto itr = cursor;

  const auto next_delim = [end_itr](auto &itr) {
    itr = std::find(itr, end_itr, '\t');
  };

  // get the read name
  const auto name_itr = itr;
  next_delim(itr);
  rec.name_len = std::distance(name_itr, itr);
  if (itr++ == end_itr)
    return false;

  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)

  // get the flag (for revcomp)
  int flag{};
  auto [flag_end, ec] = std::from_chars(itr, end_itr, flag);
  if (ec != std::errc{})
    return false;
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-const-cast)
  itr += std::distance(itr, const_cast<samrec::pos_t>(flag_end));
  if (itr++ == end_itr)
    return false;

  // skip fields
  for (auto i = 0; i < n_fields_to_skip; ++i) {
    next_delim(itr);
    if (itr++ == end_itr)
      return false;
  }

  // get the sequence
  const auto seq_itr = itr;
  next_delim(itr);
  rec.seq_len = std::distance(seq_itr, itr);
  if (itr++ == end_itr)
    return false;

  // get the quality scores
  const auto qual_itr = itr;
  next_delim(itr);
  const auto qual_len = std::distance(qual_itr, itr);
  if (itr++ == end_itr)
    return false;

  if (*qual_itr != qual_missing_symbol && qual_len != rec.seq_len)
    return false;

  const auto rec_size = rec.name_len + 2 * rec.seq_len;
  if (std::size(rec.buffer) < rec_size)
    rec.buffer.resize(rec_size);

  auto buf_beg = std::begin(rec.buffer);
  std::copy_n(name_itr, rec.name_len, buf_beg);
  if (bam_is_rev(flag)) {
    std::reverse_copy(seq_itr, seq_itr + rec.seq_len, buf_beg + rec.name_len);
    auto buf_seq_itr = buf_beg + rec.name_len;
    std::transform(buf_seq_itr, buf_seq_itr + rec.seq_len, buf_seq_itr,
                   complem);
    std::reverse_copy(qual_itr, qual_itr + rec.seq_len,
                      buf_beg + rec.name_len + rec.seq_len);
  }
  else {
    std::copy(seq_itr, seq_itr + rec.seq_len, buf_beg + rec.name_len);
    std::copy(qual_itr, qual_itr + rec.seq_len,
              buf_beg + rec.name_len + rec.seq_len);
  }
  /// ADS: removing the code below fixes a bug, and it should be sufficient, but
  /// it's not the most robust way. This needs to be revisited.
  // if (*qual_itr == qual_missing_symbol)  // could use to eliminate copy
  //   buf_beg[rec.name_len + rec.seq_len] =
  //   static_cast<char>(qual_missing_code);

  itr = std::find(itr, end_itr, '\n');
  if (*itr++ != '\n')
    return false;

  cursor = itr;

  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

  return true;
}
