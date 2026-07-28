// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_SAMREC_HPP_
#define SRC_SAMREC_HPP_

#include <cstdint>
#include <format>
#include <iterator>
#include <string>
#include <utility>
#include <vector>

#ifdef bam_is_rev
#undef bam_is_rev
#endif

#ifdef BAM_FREVERSE
#undef BAM_FREVERSE
#endif

class samrec {
public:
  using pos_t = char *;

private:
  static constexpr auto qual_missing_symbol = '*';         // from SAMv1.pdf
  static constexpr std::uint8_t qual_missing_code = 0xff;  // from sam.c

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
  // clang-format on

  // clang-format off
  samrec() = default;
  ~samrec() = default;
  samrec(const samrec &) = delete;
  auto operator=(const samrec &) -> samrec & = delete;
  auto operator=(samrec &&) noexcept -> samrec & = delete;
  samrec(samrec &&) noexcept = delete;
  // clang-format on

  [[nodiscard]] auto
  to_string() const {
    return std::format(
      "@{}\n{}\n+\n{}",
      std::string(std::data(buffer), std::data(buffer) + name_len),
      std::string(std::data(buffer) + name_len,
                  std::data(buffer) + name_len + seq_len),
      std::string(std::data(buffer) + name_len + seq_len,
                  std::data(buffer) + name_len + seq_len + seq_len));
  }

  operator bool() const { return name_len != 0; }

  [[nodiscard]] static auto
  get_next(samrec::pos_t &cursor, const samrec::pos_t end_itr,
           samrec &rec) -> bool;

  [[nodiscard]] static auto
  find_end_pos(pos_t itr, const pos_t end) -> samrec::pos_t;
};

// NOLINTBEGIN(cppcoreguidelines-pro-bounds-pointer-arithmetic)
[[nodiscard]] inline constexpr auto
get_name(const samrec &rec) {
  return std::data(rec.buffer);
}

[[nodiscard]] inline constexpr auto
get_name_end(const samrec &rec) {
  return get_name(rec) + rec.name_len;
}

[[nodiscard]] inline constexpr auto
get_seq(const samrec &rec) {
  return get_name_end(rec);
}

[[nodiscard]] inline constexpr auto
get_seq_end(const samrec &rec) {
  return get_seq(rec) + rec.seq_len;
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
  return get_qual(rec) + rec.seq_len;
}

[[nodiscard]] inline constexpr auto
get_qual_size(const samrec &rec) {
  return get_seq_size(rec);
}
// NOLINTEND(cppcoreguidelines-pro-bounds-pointer-arithmetic)

using sam_chunks_t = std::vector<std::pair<samrec::pos_t, samrec::pos_t>>;

struct sam_task_t {
  samrec::pos_t beg{};
  samrec::pos_t end{};
};

[[nodiscard]] inline auto
get_next(auto &itr, const auto end, samrec &rec) -> bool {
  return samrec::get_next(itr, end, rec);
}

#endif  // SRC_SAMREC_HPP_
