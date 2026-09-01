// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bgzf_reader.hpp"
#include "bgzf_block.hpp"

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

[[nodiscard]] static inline constexpr auto
get_unaligned_le32(const auto p) -> std::int32_t {
  std::int32_t value{};
  std::memcpy(std::addressof(value), p, sizeof(std::int32_t));
  if constexpr (std::endian::native == std::endian::big)
    return std::byteswap(value);
  else
    return value;
}

[[nodiscard]] static inline constexpr auto
get_isize(const auto data, const auto data_size) {
  static constexpr decltype(data_size) isize_size = 4;
  assert(data_size > isize_size);
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  const auto data_isize = data + data_size - isize_size;
  return data_size < isize_size ? 0 : get_unaligned_le32(data_isize);
}

auto
assign(gzip_header &hdr, const auto data) -> void {
  // ADS: data from the file takes 18 bytes, the fields of the struct take 18
  // bytes, but the struct occupies 20 due to uint32_t members
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  std::memcpy(reinterpret_cast<std::uint8_t *>(&hdr), std::to_address(data),
              gzip_header_size);
}

bgzf_reader::bgzf_reader(const std::string &filename,
                         const std::int64_t buf_size) :
  fp(std::fopen(std::data(filename), "r"), &std::fclose),     //
  filesize{std::filesystem::file_size(filename)},             //
  inbuf(std::make_unique_for_overwrite<char[]>(inbuf_size)),  // NOLINT
  outbuf(std::make_unique_for_overwrite<char[]>(buf_size)),   // NOLINT
  next_in{inbuf.get()},                                       //
  end_in{inbuf.get()},                                        //
  next_out{outbuf.get()},                                     //
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  end_out{outbuf.get() + buf_size}  //
{}

[[nodiscard]] auto
bgzf_reader::read_data() -> bool {
  if (at_eof())
    return false;
  const auto unused_in = std::distance(next_in, end_in);
  std::memcpy(inbuf.get(), next_in, unused_in);
  next_in = inbuf.get();
  // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  end_in = next_in + unused_in;
  const auto avail_in = inbuf_size - unused_in;
  const auto n_bytes = std::fread(end_in, 1, avail_in, fp.get());
  if (std::ferror(fp.get()))
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed reading input");
  // will usually be end of inbuf
  end_in += n_bytes;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return n_bytes > 0;
}

[[nodiscard]] inline constexpr auto
get_gzip_body_size(const gzip_header &gh) -> std::uint32_t {
  return (static_cast<std::uint32_t>(gh.size) + 1) - gzip_header_size;
}

[[nodiscard]] auto
bgzf_reader::get_decomp_task(char *out_itr) -> bgzf_block_t {
  if (std::distance(next_out, end_out) < max_bgzf_block_size)
    return bgzf_block_t{};
  if (std::distance(next_in, end_in) < gzip_header_size && !read_data())
    return bgzf_block_t{};
  assign(gh, next_in);
  assert(gh.check_magic());
  // NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  next_in += gzip_header_size;
  const auto body_size = get_gzip_body_size(gh);
  if (std::distance(next_in, end_in) < body_size && !read_data())
    return bgzf_block_t{};
  std::memcpy(next_out, next_in, body_size);
  bgzf_block_t task(get_isize(next_in, body_size), out_itr, next_out);
  next_in += body_size;
  next_out += max_bgzf_block_size;
  // NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)
  return task;
}
