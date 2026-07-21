// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bgzf_reader.hpp"
#include "bgzf_block.hpp"  // for bgzf_block_t, max_bgzf_block_size

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <string>
#include <system_error>

[[nodiscard]] static constexpr auto
get_unaligned_le32(const std::uint8_t *p) -> std::uint32_t {
  return (static_cast<std::uint32_t>(p[3]) << 24) |
         (static_cast<std::uint32_t>(p[2]) << 16) |
         (static_cast<std::uint32_t>(p[1]) << 8) |
         (static_cast<std::uint32_t>(p[0]) << 0);
}

[[nodiscard]] static constexpr auto
get_isize(const auto *data, const auto data_size) {
  static constexpr decltype(data_size) isize_size = 4;
  assert(data_size > isize_size);
  const auto u_data = reinterpret_cast<const std::uint8_t *>(data);
  return data_size < isize_size
           ? 0
           : get_unaligned_le32(u_data + data_size - isize_size);
}

auto
assign(gz_hdr &gz, auto *data) -> void {
  std::memcpy(reinterpret_cast<std::uint8_t *>(&gz), data, gzip_header_size);
}

bgzf_reader::bgzf_reader(const std::string &filename,
                         const std::int32_t buf_size) :
  fp(std::fopen(std::data(filename), "r"), &std::fclose), buf_size{buf_size},
  inbuf(buf_size), outbuf(buf_size), next_in{std::data(inbuf)},
  end_in{std::data(inbuf)}, next_out{std::data(outbuf)},
  end_out{std::data(outbuf) + buf_size} {}

[[nodiscard]] auto
bgzf_reader::read_data() -> bool {
  if (std::feof(fp.get()))
    return false;
  const auto unused_bytes = std::distance(next_in, end_in);
  std::copy(next_in, end_in, std::data(inbuf));
  next_in = std::data(inbuf);
  end_in = next_in + unused_bytes;
  const std::int32_t avail_in =
    buf_size - std::distance(std::data(inbuf), end_in);
  const auto r = std::fread(end_in, 1, avail_in, fp.get());
  if (std::ferror(fp.get()))
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed reading input");
  end_in += r;  // will usually be end of inbuf
  return r > 0;
}

[[nodiscard]] auto
bgzf_reader::get_decomp_task() -> bgzf_block_t {
  // ADS: need to ensure each constructed bgzf_block_t has an out_itr set and
  // not to nullptr
  if (std::distance(next_out, end_out) < max_bgzf_block_size)
    return bgzf_block_t{};
  if (std::distance(next_in, end_in) < gzip_header_size && !read_data())
    return bgzf_block_t{};
  assign(gh, next_in);
  assert(gh.check_magic());
  next_in += gzip_header_size;
  const auto body_size =
    (static_cast<std::uint32_t>(gh.size) + 1) - gzip_header_size;
  if (std::distance(next_in, end_in) < body_size && !read_data())
    return bgzf_block_t{};
  const auto uncompressed_size = get_isize(next_in, body_size);
  std::memcpy(next_out, next_in, body_size);
  bgzf_block_t task(uncompressed_size, nullptr, next_out);
  next_in += body_size;
  next_out += max_bgzf_block_size;
  return task;
}
