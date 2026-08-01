// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_SAM_FILE_HPP_
#define SRC_SAM_FILE_HPP_

#include "samrec.hpp"

#include <atomic>
#include <cstdint>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

struct task_queue;

class sam_file {
  static constexpr auto min_buf_size = 64 * 1024;
  std::vector<char> buffer;
  std::int64_t last{};
  std::int64_t cursor{};
  std::unique_ptr<std::FILE, int (*)(std::FILE *)> in;

public:
  sam_file(const std::string &filename, const std::int64_t buf_size);
  operator bool() const { return cursor < last || !std::feof(in.get()); }

  // clang-format off
  sam_file(const sam_file &) = delete;
  auto operator=(const sam_file &) -> sam_file & = delete;
  auto operator=(sam_file &&) noexcept -> sam_file & = delete;
  sam_file(sam_file &&) noexcept = default;
  ~sam_file() = default;
  // clang-format on

  auto
  make_tasks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

  auto
  reset() -> void {
    buffer.clear();
    buffer.shrink_to_fit();
  }

private:
  auto
  get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

  auto
  shift_output_buffer() -> void;

  [[nodiscard]] auto
  skip_header() -> bool;
};

inline auto
make_tasks(sam_file &reads_file,          //
           const std::int64_t n_threads,  //
           const std::int32_t file_id,    //
           task_queue &tq,                //
           std::atomic_int32_t &n_tasks   //
           ) -> void {
  static constexpr auto n_chunks_per_thread = 8;
  const auto n_chunks = n_chunks_per_thread * n_threads;
  reads_file.make_tasks(n_chunks, file_id, tq, n_tasks);
}

inline auto
reset(sam_file &reads_file) -> void {
  reads_file.reset();
}

#endif  // SRC_SAM_FILE_HPP_
