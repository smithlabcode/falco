// SPDX-License-Identifier: MIT; (c) 2026 Andrew D Smith

#include "bgzf_block.hpp"

#include <libdeflate.h>

#include <cassert>
#include <cstdlib>

auto
bgzf_block_t::decompress() -> void {
  // ADS: how much redundant allocation happens here?
  auto decompressor = libdeflate_alloc_decompressor();
  std::size_t inflated_size{};
  [[maybe_unused]]
  const libdeflate_result result = libdeflate_deflate_decompress(  //
    decompressor,                                                  //
    in_itr,                                                        //
    max_bgzf_block_size,                                           //
    out_itr,                                                       //
    size,                                                          //
    &inflated_size);
  assert(result == LIBDEFLATE_SUCCESS &&
         inflated_size == static_cast<std::size_t>(size));
  libdeflate_free_decompressor(decompressor);
}
