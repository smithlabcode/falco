/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SRC_FILE_INFO_HPP_
#define SRC_FILE_INFO_HPP_

#include "falco_file_format.hpp"

#include "nlohmann/json.hpp"

#include <cstdint>
#include <string>

struct file_info {
  std::string name;
  falco::file_format format{};
  std::string description;
  std::int64_t size{};
  std::uint64_t n_reads_est{};
  std::uint64_t read_len_est{};
  falco::encoding encoding{};
  bool has_tiles{};
  std::uint32_t tile_id_position{};

  [[nodiscard]] auto
  string() const -> std::string {
    static constexpr auto n_indent = 4;
    return nlohmann::json(*this).dump(n_indent);
  }

  auto
  set_encoding(const falco::encoding &e) {
    encoding = e;
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(file_info, name, format, description, size,
                                 n_reads_est, read_len_est, encoding, has_tiles,
                                 tile_id_position);
};

#endif  // SRC_FILE_INFO_HPP_
