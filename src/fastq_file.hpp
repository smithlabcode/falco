// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FASTQ_FILE_HPP_
#define SRC_FASTQ_FILE_HPP_

#ifdef USE_MMAP_FOR_PLAIN_FASTQ
#include <unistd.h>
#endif

#include <atomic>
#include <cstdint>
#include <string>
#include <tuple>
#include <variant>
#ifndef USE_MMAP_FOR_PLAIN_FASTQ
#include <compare>
#include <cstdio>
#include <iterator>
#include <memory>  // IWYU pragma: keep
#include <vector>
#else  // NOT USE_MMAP_FOR_PLAIN_FASTQ
#include <span>
#endif

struct task_queue;

#ifndef USE_MMAP_FOR_PLAIN_FASTQ

class fastq_file {
  static constexpr auto min_buf_size = 256 * 1024;
  std::vector<char> buffer;
  std::vector<char>::iterator cursor;
  std::vector<char>::iterator last;
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> in;

public:
  fastq_file(const std::string &filename, const std::int64_t buf_size);
  operator bool() const { return cursor < last || !std::feof(in.get()); }

  // clang-format off
  fastq_file(const fastq_file &) = delete;
  auto operator=(const fastq_file &) -> fastq_file & = delete;
  auto operator=(fastq_file &&) noexcept -> fastq_file & = delete;
  fastq_file(fastq_file &&) noexcept = default;
  ~fastq_file() = default;
  // clang-format on

  friend auto
  reset(fastq_file &reads_file) -> void;

  friend auto
  make_tasks(fastq_file &reads_file,
             const std::int64_t n_threads,
             const std::int32_t file_id,
             task_queue &tq,
             std::atomic_int32_t &n_tasks) -> void;

private:
  auto
  reset() -> void {
    buffer.clear();
    buffer.shrink_to_fit();
    cursor = std::begin(buffer);
    last = cursor;
  }

  auto
  get_chunks(const std::int64_t n_chunks,
             const std::int32_t file_id,
             task_queue &tq,
             std::atomic_int32_t &n_tasks) -> void;

  auto
  make_tasks(const std::int64_t n_chunks,
             const std::int32_t file_id,
             task_queue &tq,
             std::atomic_int32_t &n_tasks) -> void;

  auto
  shift_output_buffer() -> void;

  auto
  load_next() -> void;

  [[nodiscard]] auto
  skip_header() -> bool;
};

#else  // USE_MMAP_FOR_PLAIN_FASTQ

struct fastq_file {
  static constexpr auto min_buf_size = 256 * 1024;
  std::int64_t target_length{};
  std::int64_t filesize{};
  char *mmap_data{};
  std::int64_t length{};
  std::int64_t start_in_file{};
  std::int64_t stop_in_file{};
  int fd{};
  std::span<char> buffer{};
  std::span<char>::iterator last{};

  fastq_file(const std::string &filename, const std::int64_t target_length);

  // clang-format off
  fastq_file(const fastq_file &) = delete;
  auto operator=(const fastq_file &) -> fastq_file & = delete;
  auto operator=(fastq_file &&) noexcept -> fastq_file & = delete;
  // clang-format on

  fastq_file(fastq_file &&src) noexcept :
    target_length{src.target_length},  //
    filesize{src.filesize},            //
    mmap_data{src.mmap_data},          //
    length{src.length},                //
    start_in_file{src.start_in_file},  //
    stop_in_file{src.stop_in_file},    //
    fd{dup(src.fd)},                   // <- LOOK
    buffer{src.buffer},                //
    last{src.last}                     //
  {}

  ~fastq_file() {
    reset();
    close(fd);  // will always have been opened using a filename
  }

  operator bool() const { return stop_in_file != filesize; }

  friend auto
  reset(fastq_file &reads_file) -> void;

  friend auto
  make_tasks(fastq_file &reads_file,
             const std::int64_t n_threads,
             const std::int32_t file_id,
             task_queue &tq,
             std::atomic_int32_t &n_tasks) -> void;

private:
  auto
  reset() -> void;

  auto
  load_next() -> void;

  auto
  get_chunks(std::int64_t n_chunks,
             const std::int32_t file_id,
             task_queue &tq,
             std::atomic_int32_t &n_tasks) -> void;

  auto
  make_tasks(const std::int64_t n_threads,
             const std::int32_t file_id,
             task_queue &tq,
             std::atomic_int32_t &n_tasks) -> void;
};

#endif  // USE_MMAP_FOR_PLAIN_FASTQ

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

inline auto
make_tasks(fastq_file &reads_file,
           const std::int64_t n_threads,
           const std::int32_t file_id,
           task_queue &tq,
           std::atomic_int32_t &n_tasks) -> void {
  n_tasks = 1;  // for current task, which makes tasks
  reads_file.make_tasks(n_threads, file_id, tq, n_tasks);
}

inline auto
reset(fastq_file &reads_file) -> void {
  reads_file.reset();
}

#endif  // SRC_FASTQ_FILE_HPP_
