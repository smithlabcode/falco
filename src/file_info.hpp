// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FILE_INFO_HPP_
#define SRC_FILE_INFO_HPP_

#include "nlohmann/json.hpp"

#include <cstdint>
#include <string>

namespace falco {
enum class encoding : std::uint8_t;
enum class file_format : std::uint8_t;
}  // namespace falco

struct file_info {
  std::string name;
  falco::file_format format{};
  std::string description;
  std::int64_t size{};
  std::uint64_t n_reads_est{};
  std::uint64_t read_len_est{};
  bool has_quals{};
  falco::encoding encoding{};
  bool has_tiles{};
  std::uint32_t tile_id_position{};

  [[nodiscard]] auto
  to_string() const -> std::string {
    static constexpr auto n_indent = 4;
    return nlohmann::json(*this).dump(n_indent);
  }

  [[nodiscard]] auto
  to_line() const -> std::string;

  auto
  set_encoding(const falco::encoding &e) {
    encoding = e;
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(file_info,
                                 name,
                                 format,
                                 description,
                                 size,
                                 n_reads_est,
                                 read_len_est,
                                 has_quals,
                                 encoding,
                                 has_tiles,
                                 tile_id_position);
};

#endif  // SRC_FILE_INFO_HPP_
