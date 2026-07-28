// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BAM_FILE_HPP_
#define SRC_BAM_FILE_HPP_

#include "bam_header.hpp"
#include "bamrec.hpp"
#include "bgzf_block.hpp"
#include "bgzf_reader.hpp"

#include <atomic>
#include <cstdint>
#include <iterator>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

struct task_queue;

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
    input_buffer(buf_size / 3 + min_buf_size),   //
    output_buffer(buf_size / 3 + min_buf_size),  //
    br(filename, buf_size / 3)                   //
  {}

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
  load_next(const std::int32_t file_id, task_queue &tq,
            std::atomic_int32_t &n_tasks) -> void;

private:
  auto
  get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

public:
  [[nodiscard]] auto
  inflate_only() const -> bool {
    return is_first_load;
  }

  auto
  make_tasks_inflate(const std::int32_t file_id, task_queue &tq,
                     std::atomic_int32_t &n_tasks) -> void;

  auto
  make_tasks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

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

inline auto
make_tasks(bam_file &reads_file,         //
           const std::int64_t n_chunks,  //
           const std::int32_t file_id,   //
           task_queue &tq,               //
           std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes tasks
  if (!reads_file.inflate_only())
    reads_file.make_tasks(n_chunks, file_id, tq, n_tasks);
  reads_file.load_next(file_id, tq, n_tasks);
  reads_file.read_data();
}

#endif  // SRC_BAM_FILE_HPP_
