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

#ifndef SRC_RUN_MODE_HPP_
#define SRC_RUN_MODE_HPP_

#include "boost/boost_unordered.hpp"

#include <string>
#include <vector>

struct run_mode {
  bool do_tiles{};
  bool do_kmers{};
  bool do_dups{};
  bool do_groups{};

  auto
  assign(const boost::unordered_flat_map<std::string, bool> &modes) -> void;

  [[nodiscard]] auto
  string() const -> std::string {
    return {};
  }

  static std::vector<std::string> labels;
};

[[nodiscard]] constexpr inline auto
do_tiles(const run_mode &rm) {
  return rm.do_tiles;
}

[[nodiscard]] constexpr inline auto
do_dups(const run_mode &rm) {
  return rm.do_dups;
}

[[nodiscard]] constexpr inline auto
do_kmers(const run_mode &rm) {
  return rm.do_kmers;
}

[[nodiscard]] constexpr inline auto
do_groups(const run_mode &rm) {
  return rm.do_groups;
}

#endif  // SRC_RUN_MODE_HPP_
