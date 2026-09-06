// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FASTQ_STDIN_HPP_
#define SRC_FASTQ_STDIN_HPP_

#include <atomic>
#include <compare>
#include <cstdint>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>

struct task_queue;

struct fastq_stdin {
  static constexpr auto min_buf_size = 64 * 1024;
  std::vector<char> buffer;
  std::vector<char>::iterator cursor;
  std::vector<char>::iterator last;
  bool hit_eof{};

  explicit fastq_stdin(const std::int64_t buf_size);
  operator bool() const { return cursor < last || !hit_eof; }

  // clang-format off
  fastq_stdin(const fastq_stdin &) = delete;
  auto operator=(const fastq_stdin &) -> fastq_stdin & = delete;
  auto operator=(fastq_stdin &&) noexcept -> fastq_stdin & = delete;
  fastq_stdin(fastq_stdin &&) noexcept = default;
  ~fastq_stdin() = default;
  // clang-format on

  auto
  reset() -> void {
    buffer.clear();
    buffer.shrink_to_fit();
    cursor = std::begin(buffer);
    last = std::begin(buffer);
  }

  auto
  make_tasks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

private:
  auto
  get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

  auto
  shift_output_buffer() -> void;

  auto
  load_next() -> void;
};

[[nodiscard]] auto
estimate_n_reads_fastq_stdin(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

inline auto
make_tasks(fastq_stdin &reads_file,       //
           const std::int64_t n_threads,  //
           const std::int32_t file_id,    //
           task_queue &tq,                //
           std::atomic_int32_t &n_tasks) -> void {
  static constexpr auto n_chunks_per_thread = 8;
  const auto n_chunks = n_chunks_per_thread * n_threads;
  reads_file.make_tasks(n_chunks, file_id, tq, n_tasks);
}

inline auto
reset(fastq_stdin &reads_file) -> void {
  reads_file.reset();
}

#endif  // SRC_FASTQ_STDIN_HPP_
