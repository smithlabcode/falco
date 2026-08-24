// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bamrec.hpp"
#include "quality_score.hpp"

#include <algorithm>
#include <span>
#include <string>

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
