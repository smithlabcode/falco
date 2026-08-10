// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FASTQ_BGZF_FILE_HPP_
#define SRC_FASTQ_BGZF_FILE_HPP_

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

class fastq_bgzf_file {
  static constexpr std::int64_t min_buf_size = 64 * 1024;
  std::vector<char> input_buffer;
  std::vector<char> output_buffer;
  std::int64_t input_last{};
  std::int64_t output_last{};
  std::int64_t input_cursor{};
  std::int64_t output_cursor{};
  bgzf_reader br;
  bool is_first_load{true};
  bool had_last_chunks{false};

  [[nodiscard]] static auto
  get_reader_buffer_size(const std::int64_t buf_size) -> std::int64_t {
    return buf_size / 3;  // 1/3 here seems to work well...
  }

  [[nodiscard]] static auto
  get_input_buffer_size(const std::int64_t buf_size) -> std::int64_t {
    // split buffer capacity between input and output buffers
    return (buf_size - get_reader_buffer_size(buf_size)) / 2;
  }

  [[nodiscard]] static auto
  get_output_buffer_size(const std::int64_t buf_size) -> std::int64_t {
    // split buffer capacity between input and output buffers
    return (buf_size - get_reader_buffer_size(buf_size)) / 2;
  }

public:
  fastq_bgzf_file(const std::string &filename, const std::int64_t buf_size) :
    // Three different buffers, 1/4 goes to the bgzf_reader
    input_buffer(get_input_buffer_size(buf_size) + min_buf_size),    //
    output_buffer(get_output_buffer_size(buf_size) + min_buf_size),  //
    br(filename, get_reader_buffer_size(buf_size))                   //
  {}

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

  friend auto
  reset(fastq_bgzf_file &reads_file) -> void;

  friend auto
  make_tasks(fastq_bgzf_file &reads_file,   //
             const std::int64_t n_threads,  //
             const std::int32_t file_id,    //
             task_queue &tq,                //
             std::atomic_int32_t &n_tasks) -> void;

private:
  auto
  reset() {
    input_buffer.clear();
    input_buffer.shrink_to_fit();
    output_buffer.clear();
    output_buffer.shrink_to_fit();
    br.release();
  }

  auto
  get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

  auto
  load_next(const std::int32_t file_id, task_queue &tq,
            std::atomic_int32_t &n_tasks) -> void;

  auto
  make_tasks_inflate(const std::int32_t file_id, task_queue &tq,
                     std::atomic_int32_t &n_tasks) -> void;

  auto
  make_tasks(const std::int64_t n_threads, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

  [[nodiscard]] auto
  inflate_only() const -> bool {
    return is_first_load;
  }

  auto
  read_data() -> void {
    [[maybe_unused]] const auto r = br.read_data();
  }

  [[nodiscard]] auto
  has_in() const -> bool {
    return input_last + max_bgzf_block_size < std::ssize(input_buffer);
  }
};

[[nodiscard]] auto
estimate_n_reads_fastq_bgzf(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

inline auto
make_tasks(fastq_bgzf_file &reads_file,   //
           const std::int64_t n_threads,  //
           const std::int32_t file_id,    //
           task_queue &tq,                //
           std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // +1 so not to decrement to zero until current thread done
  if (!reads_file.inflate_only())
    reads_file.make_tasks(n_threads, file_id, tq, n_tasks);
  reads_file.load_next(file_id, tq, n_tasks);
  reads_file.read_data();
}

inline auto
reset(fastq_bgzf_file &reads_file) -> void {
  reads_file.reset();
}

#endif  // SRC_FASTQ_BGZF_FILE_HPP_
