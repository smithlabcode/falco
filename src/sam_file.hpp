// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_SAM_FILE_HPP_
#define SRC_SAM_FILE_HPP_

#include "bamrec.hpp"

#include <atomic>
#include <cstdint>
#include <iterator>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

struct task_queue;

class sam_file {
public:
  using rec_t = bamrec;
  static constexpr auto min_buf_size = 64 * 1024;
  std::vector<char> buffer;
  std::int64_t last{};
  std::int64_t cursor{};

  sam_file([[maybe_unused]] const std::string &filename,
           const std::int64_t buf_size) : buffer(buf_size + min_buf_size) {
    throw std::runtime_error("Not implemented. Workaround: make BAM files");
  }

  // clang-format off
  // delete copy and assignment
  sam_file(const sam_file &) = delete;
  auto operator=(const sam_file &) -> sam_file & = delete;
  auto operator=(sam_file &&) noexcept -> sam_file & = delete;
  // default move for emplace
  sam_file(sam_file &&) noexcept = default;
  ~sam_file() = default;
  // clang-format on

  [[nodiscard]] operator bool() const { return false; }

  auto
  load_next(const std::int32_t file_id, task_queue &tq,
            std::atomic_int32_t &n_tasks) -> void;

private:
  auto
  get_chunks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;

public:
  auto
  make_tasks(const std::int64_t n_chunks, const std::int32_t file_id,
             task_queue &tq, std::atomic_int32_t &n_tasks) -> void;
};

[[nodiscard]] inline constexpr auto
is_active(const sam_file &reads_file) -> bool {
  return static_cast<bool>(reads_file);
}

inline auto
make_tasks([[maybe_unused]] sam_file &reads_file,         //
           [[maybe_unused]] const std::int64_t n_chunks,  //
           [[maybe_unused]] const std::int32_t file_id,   //
           [[maybe_unused]] task_queue &tq,               //
           [[maybe_unused]] std::atomic_int32_t &n_tasks  //
           ) -> void {
  n_tasks = 1;  // for current task, which makes tasks
  // work goes here
}

#endif  // SRC_SAM_FILE_HPP_
