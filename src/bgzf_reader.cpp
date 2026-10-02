// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bgzf_reader.hpp"

#include "bgzf_block.hpp"

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>
#include <utility>

bgzf_reader::bgzf_reader(const std::string &filename,
                         const std::int64_t outbuf_size) :
  fp(std::fopen(std::data(filename), "r"), &std::fclose),  //
  filesize{std::filesystem::file_size(filename)},          //
  inbuf(inbuf_size),                                       //
  outbuf(outbuf_size),                                     //
  next_in_itr{std::begin(inbuf)},                          //
  end_in_itr{std::begin(inbuf)},                           //
  next_out_itr{std::begin(outbuf)},                        //
  end_out_itr{std::end(outbuf)}                            //
{}

[[nodiscard]] auto
bgzf_reader::read_data() -> bool {
  if (at_eof())
    return false;
  const auto unused_in = std::distance(next_in_itr, end_in_itr);
  std::copy_n(next_in_itr, unused_in, std::begin(inbuf));
  next_in_itr = std::begin(inbuf);
  end_in_itr = next_in_itr + unused_in;
  const auto avail_in = inbuf_size - unused_in;
  const auto n_bytes =
    std::fread(std::to_address(end_in_itr), 1, avail_in, fp.get());
  if (std::ferror(fp.get()))
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed reading input");
  end_in_itr += n_bytes;  // will usually be end of inbuf
  return n_bytes > 0;
}

namespace {

static constexpr auto gzip_header_size = 18;

struct gzip_header {
  static constexpr auto magic1 = 0x1F;
  static constexpr auto magic2 = 0x8B;
  static constexpr auto size_index = 16;

  std::array<std::uint8_t, gzip_header_size> data{};

  /// internal structure
  // std::uint8_t id1{};       // 0
  // std::uint8_t id2{};       // 1
  // std::uint8_t cm_eight{};  // 2
  // std::uint8_t flg{};       // 3
  // std::uint32_t mtime{};    // 7 [4]
  // std::uint8_t xfl{};       // 8
  // std::uint8_t os{};        // 9
  // std::uint16_t xlen{};     // 11 [2]
  // char b{};                 // 12
  // char c{};                 // 13
  // std::uint16_t two{};      // 14 [2]
  // std::uint16_t size{};     // 16 [2]

  [[nodiscard]] auto
  get_size() const -> std::uint16_t {
    std::uint16_t s{};
    std::memcpy(std::addressof(s), std::addressof(data[size_index]),
                sizeof(std::uint16_t));
    return s;
  }

  [[nodiscard]] auto
  check_magic() const -> bool {
    return data[0] == magic1 && data[1] == magic2;  // ADS: check b and c also?
  }
};

[[nodiscard]] static inline constexpr auto
get_unaligned_le32(const auto p) -> std::int32_t {
  std::int32_t value{};
  std::memcpy(std::addressof(value), std::to_address(p), sizeof(std::int32_t));
  if constexpr (std::endian::native == std::endian::big)
    return std::byteswap(value);
  else
    return value;
}

[[nodiscard]] static inline constexpr auto
get_isize(const auto data, const auto data_size) {
  static constexpr decltype(data_size) isize_size = 4;
  assert(data_size > isize_size);
  const auto data_isize = data + data_size - isize_size;
  return data_size < isize_size ? 0 : get_unaligned_le32(data_isize);
}

static inline auto
assign_gzip_hdr(gzip_header &hdr,  // cppcheck-suppress constParameterReference
                const auto data) -> void {
  std::memcpy(std::data(hdr.data), std::to_address(data), gzip_header_size);
}

[[nodiscard]] static inline constexpr auto
get_gzip_body_size(const gzip_header &gh) -> std::uint32_t {
  return (static_cast<std::uint32_t>(gh.get_size()) + 1) - gzip_header_size;
}

}  // namespace

[[nodiscard]] auto
bgzf_reader::get_decomp_task(bgzf_reader::iterator out_itr) -> bgzf_block_t {
  if (std::distance(next_out_itr, end_out_itr) < max_bgzf_block_size)
    return bgzf_block_t{};
  if (std::distance(next_in_itr, end_in_itr) < gzip_header_size && !read_data())
    return bgzf_block_t{};
  gzip_header gh;
  assign_gzip_hdr(gh, next_in_itr);
  assert(gh.check_magic());
  next_in_itr += gzip_header_size;
  const auto body_size = get_gzip_body_size(gh);
  if (std::distance(next_in_itr, end_in_itr) < body_size && !read_data())
    return bgzf_block_t{};
  std::copy_n(next_in_itr, body_size, next_out_itr);
  const auto isize = get_isize(next_in_itr, body_size);
  assert(isize >= 0);
  bgzf_block_t task(isize, out_itr, next_out_itr);
  next_in_itr += body_size;
  next_out_itr += max_bgzf_block_size;
  return task;
}
