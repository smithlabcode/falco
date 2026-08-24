// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FASTQ_FILE_HPP_
#define SRC_FASTQ_FILE_HPP_

#include <fcntl.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <cstdint>
#include <filesystem>
#include <iterator>
#include <string>
#include <system_error>
#include <tuple>
#include <variant>

struct task_queue;

struct fastq_file {
  static constexpr auto min_buf_size = 65536;
  std::int64_t target_buffer_size{};
  std::int64_t filesize{};
  char *buffer{};
  std::int64_t buffer_size{};
  std::int64_t start_in_file{};
  std::int64_t stop_in_file{};
  std::int64_t cursor{};
  int fd{};

  fastq_file(const std::string &filename,
             const std::int64_t target_buffer_size) :
    target_buffer_size{
      std::max(target_buffer_size, static_cast<std::int64_t>(min_buf_size))},
    filesize{static_cast<std::int64_t>(std::filesystem::file_size(filename))},
    fd{open(std::data(filename), O_RDONLY, 0)} {
    if (fd < 0)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
  }

  // clang-format off
  fastq_file(const fastq_file &) = delete;
  auto operator=(const fastq_file &) -> fastq_file & = delete;
  auto operator=(fastq_file &&) noexcept -> fastq_file & = delete;
  // clang-format on

  fastq_file(fastq_file &&src) noexcept :
    target_buffer_size{src.target_buffer_size},  //
    filesize{src.filesize},                      //
    buffer{src.buffer},                          //
    buffer_size{src.buffer_size},                //
    start_in_file{src.start_in_file},            //
    stop_in_file{src.stop_in_file},              //
    cursor{src.cursor},                          //
    fd{dup(src.fd)}                              // <- LOOK
  {}

  auto
  reset() -> void;

  ~fastq_file() {
    reset();
    close(fd);  // will always have been opened using a filename
  }

  operator bool() const { return stop_in_file != filesize; }

  auto
  load_next() -> void;

  auto
  get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;
};

[[nodiscard]] auto
estimate_n_reads_fastq(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

inline auto
make_tasks(fastq_file &reads_file,        //
           const std::int64_t n_threads,  //
           const std::int32_t file_id,    //
           task_queue &tq,                //
           std::atomic_int32_t &n_tasks) -> void {
  static constexpr auto n_chunks_per_thread = 8;
  const auto n_chunks = n_chunks_per_thread * n_threads;
  n_tasks = 1;  // for current task, which makes tasks
  reads_file.load_next();
  reads_file.get_chunks(n_chunks, file_id, tq, n_tasks);
}

inline auto
reset(fastq_file &reads_file) -> void {
  reads_file.reset();
}

#endif  // SRC_FASTQ_FILE_HPP_
