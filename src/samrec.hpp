// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_SAMREC_HPP_
#define SRC_SAMREC_HPP_

#include <cstdint>
#include <iterator>
#include <string>
#include <vector>

class samrec {
public:
  using pos_t = std::vector<char>::const_iterator;

private:
  static constexpr auto qual_missing_symbol = '*';  // from SAMv1.pdf
  static constexpr char qual_missing_code = -1;     // 0xff from sam.c

  std::vector<char> buffer;
  std::uint32_t name_len{};
  std::uint32_t seq_len{};

public:
  // clang-format off
  friend constexpr auto get_name(const samrec &);
  friend constexpr auto get_name_end(const samrec &);
  friend constexpr auto get_seq(const samrec &);
  friend constexpr auto get_seq_end(const samrec &);
  friend constexpr auto get_seq_size(const samrec &);
  friend constexpr auto get_qual(const samrec &);
  friend constexpr auto get_qual_end(const samrec &);
  friend constexpr auto get_qual_size(const samrec &);

  samrec() = default;
  ~samrec() = default;
  samrec(const samrec &) = delete;
  auto operator=(const samrec &) -> samrec & = delete;
  auto operator=(samrec &&) noexcept -> samrec & = delete;
  samrec(samrec &&) noexcept = delete;
  // clang-format on

  [[nodiscard]] auto
  to_string() const -> std::string;

  operator bool() const { return name_len != 0; }

  [[nodiscard]] static auto
  get_next(pos_t &cursor, const pos_t end_itr, samrec &rec) -> bool;

  [[nodiscard]] static auto
  find_end_pos(pos_t itr, const pos_t end) -> pos_t;
};

[[nodiscard]] inline constexpr auto
get_name(const samrec &rec) {
  // ADS: something is unhappy when this function returns an interator
  return std::data(rec.buffer);
}

[[nodiscard]] inline constexpr auto
get_name_end(const samrec &rec) {
  return get_name(rec) + rec.name_len;  // NOLINT
}

[[nodiscard]] inline constexpr auto
get_seq(const samrec &rec) {
  return get_name_end(rec);
}

[[nodiscard]] inline constexpr auto
get_seq_end(const samrec &rec) {
  return get_seq(rec) + rec.seq_len;  // NOLINT
}

[[nodiscard]] inline constexpr auto
get_seq_size(const samrec &rec) {
  return rec.seq_len;
}

[[nodiscard]] inline constexpr auto
get_qual(const samrec &rec) {
  return get_seq_end(rec);
}

[[nodiscard]] inline constexpr auto
get_qual_end(const samrec &rec) {
  return get_qual(rec) + rec.seq_len;  // NOLINT
}

[[nodiscard]] inline constexpr auto
get_qual_size(const samrec &rec) {
  return get_seq_size(rec);
}

struct sam_task_t {
  samrec::pos_t beg{};
  samrec::pos_t end{};
};

[[nodiscard]] inline auto
get_next(samrec::pos_t &itr, const samrec::pos_t end, samrec &rec) -> bool {
  return samrec::get_next(itr, end, rec);
}

#endif  // SRC_SAMREC_HPP_
