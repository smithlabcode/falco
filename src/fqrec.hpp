// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FQREC_HPP_
#define SRC_FQREC_HPP_

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <concepts>
#include <cstdint>
#include <cstdio>  // IWYU pragma: keep
#include <cstdlib>
#include <cstring>  // IWYU pragma: keep
#include <filesystem>
#include <format>  // IWYU pragma: keep
#include <iterator>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

// clang-format off
struct fqrec {
  using pos_t = char *;
  pos_t n{};  // start of "name"
  pos_t r{};  // start of "read"
  pos_t o{};  // start of "other"
  pos_t q{};  // start of "quality" scores
  pos_t e{};  // end of the record
  [[nodiscard]] auto
  size() const { return static_cast<std::int32_t>(std::distance(r, o)) - 1; }
  [[nodiscard]] operator bool() const { return n != nullptr; }
  [[nodiscard]] auto
  string() const -> std::string { return {n, e}; }
};
[[nodiscard]] constexpr auto get_name(const fqrec &rec) { return rec.n; }
[[nodiscard]] constexpr auto get_name_end(const fqrec &rec) { return rec.r - 1; }
[[nodiscard]] constexpr auto get_seq(const fqrec &rec) { return rec.r; }
[[nodiscard]] constexpr auto get_seq_end(const fqrec &rec) { return rec.o - 1; }
[[nodiscard]] constexpr auto get_seq_size(const fqrec &rec) { return std::size(rec); }
[[nodiscard]] constexpr auto get_qual(const fqrec &rec) { return rec.q; }
[[nodiscard]] constexpr auto get_qual_end(const fqrec &rec) { return rec.e - 1; }
[[nodiscard]] constexpr auto get_qual_size(const fqrec &rec) { return std::size(rec); }
// clang-format on

struct fq_task_t {
  fqrec::pos_t beg{};
  fqrec::pos_t end{};
};

using fq_chunks_t = std::vector<std::pair<fqrec::pos_t, fqrec::pos_t>>;

[[nodiscard]] inline auto
get_next(fqrec::pos_t &cursor, const fqrec::pos_t end_itr) -> fqrec {
  // ADS: need to make sure cursor < end_itr or we will move past
  const auto n = cursor;
  auto itr = n + 1;  // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)

  const auto next_newline = [end_itr](auto &itr) {
    itr = std::find(itr, end_itr, '\n');
  };

  // clang-format off
  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto r = itr;

  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto o = itr;

  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto q = itr;

  next_newline(itr);
  if (itr++ == end_itr) return {};
  const auto e = itr;
  // clang-format on

  cursor = e;
  return {n, r, o, q, e};
}

#endif  // SRC_FQREC_HPP_
