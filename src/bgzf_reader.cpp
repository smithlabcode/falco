// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bgzf_reader.hpp"
#include "bgzf_block.hpp"  // for bgzf_block_t, max_bgzf_block_size

#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>

// NOLINTBEGIN(*-bounds-pointer-arithmetic,*-avoid-magic-numbers,*-type-reinterpret-cast)
[[nodiscard]] static inline constexpr auto
get_unaligned_le32(const std::uint8_t *p) -> std::uint32_t {
  return (static_cast<std::uint32_t>(p[3]) << 24) |
         (static_cast<std::uint32_t>(p[2]) << 16) |
         (static_cast<std::uint32_t>(p[1]) << 8) |
         (static_cast<std::uint32_t>(p[0]) << 0);
}

[[nodiscard]] static inline constexpr auto
get_isize(const auto *data, const auto data_size) {
  static constexpr decltype(data_size) isize_size = 4;
  assert(data_size > isize_size);
  const auto u_data = reinterpret_cast<const std::uint8_t *>(data);
  const auto u_data_isize = u_data + data_size - isize_size;
  return data_size < isize_size
           ? 0
           : static_cast<std::int32_t>(get_unaligned_le32(u_data_isize));
}

auto
assign(gzip_header &hdr, auto *data) -> void {
  // NOLINTNEXTLINE(*-type-reinterpret-cast)
  std::memcpy(reinterpret_cast<std::uint8_t *>(&hdr), data, gzip_header_size);
}

bgzf_reader::bgzf_reader(const std::string &filename,
                         const std::int32_t buf_size) :
  fp(std::fopen(std::data(filename), "r"), &std::fclose),     //
  filesize{std::filesystem::file_size(filename)},             //
  inbuf(std::make_unique_for_overwrite<char[]>(inbuf_size)),  // NOLINT
  outbuf(std::make_unique_for_overwrite<char[]>(buf_size)),   // NOLINT
  next_in{inbuf.get()},                                       //
  end_in{inbuf.get()},                                        //
  next_out{outbuf.get()},                                     //
  end_out{outbuf.get() + buf_size}                            //
{}

[[nodiscard]] auto
bgzf_reader::read_data() -> bool {
  if (at_eof())
    return false;
  const auto unused_in = std::distance(next_in, end_in);
  std::memcpy(inbuf.get(), next_in, unused_in);
  next_in = inbuf.get();
  end_in = next_in + unused_in;
  const auto avail_in = inbuf_size - unused_in;
  const auto r = std::fread(end_in, 1, avail_in, fp.get());
  if (std::ferror(fp.get())) {
    const auto errc = std::make_error_code(std::errc(errno));
    throw std::system_error(errc, "failed reading input");
  }
  end_in += r;  // will usually be end of inbuf
  return r > 0;
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
  next_in += gzip_header_size;
  const auto body_size = get_gzip_body_size(gh);
  if (std::distance(next_in, end_in) < body_size && !read_data())
    return bgzf_block_t{};
  std::memcpy(next_out, next_in, body_size);
  bgzf_block_t task(get_isize(next_in, body_size), out_itr, next_out);
  next_in += body_size;
  next_out += max_bgzf_block_size;
  return task;
}
// NOLINTEND(*-bounds-pointer-arithmetic,*-avoid-magic-numbers,*-type-reinterpret-cast)
