// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "samrec.hpp"
#include "falco_utils.hpp"
#include "quality_score.hpp"

#include <algorithm>
#include <cassert>
#include <format>
#include <iterator>
#include <span>
#include <stdexcept>
#include <system_error>

#ifdef BAM_FREVERSE
#undef BAM_FREVERSE
#endif

#ifdef bam_is_rev
#undef bam_is_rev
#endif

[[nodiscard]] auto
samrec::to_string() const -> std::string {
  const auto buffer_s = std::span(std::cbegin(buffer), std::cend(buffer));
  const auto name_s = buffer_s.subspan(0, name_len);
  const auto seq_s = buffer_s.subspan(name_len, seq_len);
  const auto qual_s = buffer_s.subspan(name_len + seq_len, seq_len);
  return std::format("@{}\n{}\n+\n{}",  //
                     std::string(std::cbegin(name_s), std::cend(name_s)),
                     std::string(std::cbegin(seq_s), std::cend(seq_s)),
                     std::string(std::cbegin(qual_s), std::cend(qual_s)));
}

[[nodiscard]] auto
samrec::get_next(pos_t &cursor, const pos_t end_itr, samrec &rec) -> bool {
  static constexpr auto n_fields_to_skip = 7;
  constexpr auto complem = [](const auto x) {
    return "TNGNNNCNNNNNNNNNNNNA"[x - 'A'];
  };
  static constexpr auto bam_is_rev = [](const auto flag) -> bool {
    static constexpr auto BAM_FREVERSE = 16;
    return (flag & BAM_FREVERSE) != 0;
  };

  const auto next_delim = [end_itr](auto &itr) {
    itr = std::find(itr, end_itr, '\t');
  };

  auto itr = cursor;

  // get the read name
  const auto name_itr = itr;
  next_delim(itr);
  const auto name_len = std::distance(name_itr, itr);
  if (itr++ == end_itr)
    return false;
  rec.name_len = name_len;

  // get the flag (for revcomp)
  int flag{};
  auto [flag_end, ec] = falco::from_chars(itr, end_itr, flag);
  if (ec != std::errc{})
    return false;
  itr += std::distance(itr, flag_end);
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
  const auto seq_len = std::distance(seq_itr, itr);
  if (itr++ == end_itr)
    return false;
  rec.seq_len = seq_len;

  // get the quality scores
  const auto qual_itr = itr;
  next_delim(itr);
  const auto qual_len = std::distance(qual_itr, itr);
  if (itr++ == end_itr)
    return false;

  if (qual_len != rec.seq_len)
    throw std::runtime_error("quality scores invalid for record: " +
                             std::string(name_itr, name_itr + name_len));

  const auto rec_size = rec.name_len + 2 * rec.seq_len;
  if (std::ssize(rec.buffer) < rec_size)
    rec.buffer.resize(rec_size);

  auto out_name_itr = std::begin(rec.buffer);
  auto out_seq_itr = out_name_itr + name_len;
  auto out_qual_itr = out_seq_itr + seq_len;
  std::copy_n(name_itr, name_len, out_name_itr);

  if (bam_is_rev(flag)) {
    std::reverse_copy(seq_itr, seq_itr + seq_len, out_seq_itr);
    std::transform(out_seq_itr, out_seq_itr + seq_len, out_seq_itr, complem);
    std::reverse_copy(qual_itr, qual_itr + seq_len, out_qual_itr);
  }
  else {
    std::copy(seq_itr, seq_itr + seq_len, out_seq_itr);
    std::copy(qual_itr, qual_itr + seq_len, out_qual_itr);
  }
  assert(std::all_of(out_qual_itr, out_qual_itr + seq_len, [](const auto q) {
    return q >= 0 && q <= falco::max_qual_val;
  }));

  itr = std::find(itr, end_itr, '\n');
  if (*itr++ != '\n')
    return false;

  cursor = itr;

  return true;
}
