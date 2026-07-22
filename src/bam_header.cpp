// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bam_header.hpp"

#include <stdexcept>

[[nodiscard]] auto
bam_header::update(bam_header::iterator itr,
                   const bam_header::iterator end) -> bam_header::iterator {
  static constexpr auto msg = "incorrect BAM magic identified: {} at {}";
  const auto update_u32 = [](auto &val, const auto inc, const auto the_byte) {
    val = val | ((inc & 0xFF) << (32 - 8 * the_byte));
  };
  while (itr != end && magic_bytes_remaining) {
    const auto magic_pos = magic_size - magic_bytes_remaining;
    if (magic[magic_pos] != *itr)
      throw std::runtime_error(std::format(msg, *itr, magic_pos));
    --magic_bytes_remaining;
    ++itr;
  }
  while (itr != end && l_text_bytes_remaining) {
    update_u32(l_text, *itr, l_text_bytes_remaining);
    --l_text_bytes_remaining;
    ++itr;
  }
  while (itr != end && l_text) {
    ++itr;
    --l_text;
  }
  while (itr != end && n_ref_bytes_remaining) {
    update_u32(n_ref, *itr, n_ref_bytes_remaining);
    --n_ref_bytes_remaining;
    ++itr;
  }
  while (itr != end && n_ref) {
    while (itr != end && l_name_bytes_remaining) {
      update_u32(name_bytes_remaining, *itr, l_name_bytes_remaining);
      --l_name_bytes_remaining;
      ++itr;
    }
    while (itr != end && name_bytes_remaining) {
      --name_bytes_remaining;
      ++itr;
    }
    while (itr != end && l_ref_bytes_remaining) {
      --l_ref_bytes_remaining;
      ++itr;
    }
    if (ref_complete())  // shouldn't care about itr
      decrement_ref();
  }
  return itr;
}
