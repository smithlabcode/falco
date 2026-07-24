// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BAM_FILE_HPP_
#define SRC_BAM_FILE_HPP_

#include "bam_header.hpp"
#include "bamrec.hpp"
#include "bgzf_block.hpp"
#include "bgzf_reader.hpp"
#include "falco_task.hpp"

#include <cstdint>
#include <format>
#include <iterator>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

class bam_file {
public:
  using rec_t = bamrec;
  static constexpr auto min_buf_size = 64 * 1024;
  std::vector<char> input_buffer;
  std::vector<char> output_buffer;
  std::int64_t input_last{};
  std::int64_t output_last{};
  std::int64_t output_cursor{};
  bgzf_reader br;
  bam_header bh;
  bool is_first_load{true};
  bool had_last_chunks{false};

  bam_file(const std::string &filename, const std::int64_t buf_size) :
    // ADS: diving the buffer size here because there are multiple and the user
    // presumably gives a total
    input_buffer(buf_size / 4 + min_buf_size),
    output_buffer(buf_size / 4 + min_buf_size), br(filename, buf_size / 4) {}

  // clang-format off
  // delete copy and assignment
  bam_file(const bam_file &) = delete;
  auto operator=(const bam_file &) -> bam_file & = delete;
  auto operator=(bam_file &&) noexcept -> bam_file & = delete;
  // default move for emplace
  bam_file(bam_file &&) noexcept = default;
  ~bam_file() = default;
  // clang-format on

  [[nodiscard]] operator bool() const { return !had_last_chunks; }

  auto
  shift_buffers() -> void;

  auto
  load_next() -> std::vector<bgzf_block_t>;

  [[nodiscard]] auto
  get_chunks(const std::int64_t n_chunks) -> bam_chunks_t;

  [[nodiscard]] auto
  inflate_only() const -> bool {
    return is_first_load;
  }

  [[nodiscard]] auto
  make_tasks_inflate() -> std::vector<task_t>;

  [[nodiscard]] auto
  make_tasks(const std::int64_t n_chunks) -> std::vector<task_t>;

  auto
  read_data() -> void {
    [[maybe_unused]] const auto r = br.read_data();
  }

private:
  [[nodiscard]] auto
  has_in() const -> bool {
    return input_last + max_bgzf_block_size < std::ssize(input_buffer);
  }
};

[[nodiscard]] auto
estimate_n_reads_bam(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

[[nodiscard]] inline constexpr auto
is_active(const bam_file &reads_file) -> bool {
  return static_cast<bool>(reads_file);
}

[[nodiscard]] inline auto
inflate_only(const bam_file &reads_file) -> bool {
  return reads_file.inflate_only();
}

[[nodiscard]] inline auto
make_tasks_inflate(bam_file &reads_file) -> std::vector<task_t> {
  return reads_file.make_tasks_inflate();
}

[[nodiscard]] inline auto
make_tasks(bam_file &reads_file,
           const std::int64_t n_chunks) -> std::vector<task_t> {
  return reads_file.make_tasks(n_chunks);
}

inline auto
read_data(bam_file &reads_file) -> void {
  reads_file.read_data();
}

#endif  // SRC_BAM_FILE_HPP_
