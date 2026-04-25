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

#ifndef SRC_FASTQ_RECORD_HPP_
#define SRC_FASTQ_RECORD_HPP_

#include "fastq_buffer.hpp"

#include <cstdint>
#include <format>
#include <limits>
#include <numeric>
#include <string>

struct fqrec {
  using fqrec_pos_t = std::int64_t;
  static constexpr auto sentinel = std::numeric_limits<fqrec_pos_t>::max();

  fqrec_pos_t n{};  // start of "name"
  fqrec_pos_t r{};  // start of "read"
  fqrec_pos_t o{};  // start of "other"
  fqrec_pos_t q{};  // start of "quality" scores
  fqrec_pos_t e{};  // end of the record

  [[nodiscard]] auto
  size() const -> std::uint32_t {
    return o - r - 1;
  }

  [[nodiscard]] operator bool() const {
    return e != std::numeric_limits<fqrec_pos_t>::max();
  }

  [[nodiscard]] auto
  string(const char *data) const -> std::string {
    return {data + n, data + e};
  }
};
static constexpr fqrec null_rec = {0, 0, 0, 0, fqrec::sentinel};

// clang-format off
[[nodiscard]] constexpr auto
get_name(const fastq_buffer &b, const fqrec &rec) { return b.data + rec.n; }
[[nodiscard]] constexpr auto
get_name_size(const fqrec &rec) { return rec.r - rec.n - 1; }
[[nodiscard]] constexpr auto
get_seq(const fastq_buffer &b, const fqrec &rec) { return b.data + rec.r; }
[[nodiscard]] constexpr auto
get_seq_size(const fqrec &rec) { return std::size(rec); }
[[nodiscard]] constexpr auto
get_qual(const fastq_buffer &b, const fqrec &rec) { return b.data + rec.q; }
[[nodiscard]] constexpr auto
get_qual_size(const fqrec &rec) { return std::size(rec); }
// clang-format on

[[nodiscard]] inline auto
get_next(const auto data, std::int64_t &cursor,
         const std::int64_t lim) -> fqrec {
  const auto n = cursor;
  auto itr = data + n;  // should point to '@', should not be '\n'
  const auto fq_end = data + lim;

  itr = std::find(itr + 1, fq_end, '\n');
  if (itr == fq_end)
    return null_rec;
  const auto r = std::distance(data, itr) + 1;

  itr = std::find(itr + 1, fq_end, '\n');
  if (itr == fq_end)
    return null_rec;
  const auto o = std::distance(data, itr) + 1;

  itr = std::find(itr + 1, fq_end, '\n');
  if (itr == fq_end)
    return null_rec;
  const auto q = std::distance(data, itr) + 1;

  itr = std::find(itr + 1, fq_end, '\n');
  if (itr == fq_end)
    return null_rec;
  const auto e = std::distance(data, itr) + 1;

  cursor = e;

  return {n, r, o, q, e};
}

#endif  // SRC_FASTQ_RECORD_HPP_
