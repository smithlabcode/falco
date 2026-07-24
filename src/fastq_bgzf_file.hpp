// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FASTQ_BGZF_FILE_HPP_
#define SRC_FASTQ_BGZF_FILE_HPP_

#include "bgzf_block.hpp"
#include "bgzf_reader.hpp"
#include "falco_task.hpp"
#include "fqrec.hpp"

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

class fastq_bgzf_file {
  using rec_t = fqrec;
  static constexpr auto min_buf_size = 64 * 1024;
  std::int64_t buf_size{};  // size of allocated buffer
  std::vector<char> input_buffer;
  std::vector<char> output_buffer;
  std::int64_t input_last{};
  std::int64_t output_last{};
  std::int64_t input_cursor{};
  std::int64_t output_cursor{};
  bgzf_reader br;
  bool is_first_load{true};
  bool had_last_chunks{false};

public:
  fastq_bgzf_file(const std::string &filename, const std::int64_t buf_size) :
    buf_size{buf_size / 2}, input_buffer(buf_size / 2),
    output_buffer(buf_size / 2), br(filename, buf_size / 2) {}

  // clang-format off
  // delete copy and assignment
  fastq_bgzf_file(const fastq_bgzf_file &) = delete;
  auto operator=(const fastq_bgzf_file &) -> fastq_bgzf_file & = delete;
  auto operator=(fastq_bgzf_file &&) noexcept -> fastq_bgzf_file & = delete;
  // default move for emplace
  fastq_bgzf_file(fastq_bgzf_file &&) noexcept = default;
  ~fastq_bgzf_file() = default;
  // clang-format on

  [[nodiscard]] operator bool() const { return !had_last_chunks; }

  auto
  load_next() -> std::vector<bgzf_block_t>;

  [[nodiscard]] auto
  get_chunks(const std::int64_t n_chunks) -> fq_chunks_t;

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
  auto
  shift_buffers() -> void;

  [[nodiscard]] auto
  has_in() const -> bool {
    return input_last + max_bgzf_block_size < std::ssize(input_buffer);
  }
};

[[nodiscard]] inline auto
inflate_only(fastq_bgzf_file &f) -> bool {
  return f.inflate_only();
}

[[nodiscard]] auto
get_chunks_fastq(fastq_bgzf_file &fq,
                 const std::int64_t n_chunks) -> fq_chunks_t;

[[nodiscard]] auto
estimate_n_reads_fastq_bgzf(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

[[nodiscard]] inline auto
make_tasks(fastq_bgzf_file &reads_file,
           const std::int64_t n_chunks) -> std::vector<task_t> {
  return reads_file.make_tasks(n_chunks);
}

[[nodiscard]] inline auto
make_tasks_inflate(fastq_bgzf_file &reads_file) -> std::vector<task_t> {
  return reads_file.make_tasks_inflate();
}

[[nodiscard]] inline auto
is_active(const fastq_bgzf_file &reads_file) -> bool {
  return static_cast<bool>(reads_file);
}

inline auto
read_data(fastq_bgzf_file &reads_file) -> void {
  reads_file.read_data();
}

#endif  // SRC_FASTQ_BGZF_FILE_HPP_
