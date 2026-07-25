// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_READS_FILE_HPP_
#define SRC_READS_FILE_HPP_

#include "task_queue.hpp"

#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <utility>  // IWYU pragma: keep
#include <variant>
#include <vector>

class reads_file_t {
public:
  // clang-format off
  template <typename T> explicit reads_file_t(T x) : self_(new model<T>(std::move(x))) {}
  reads_file_t(const reads_file_t &x) = delete;
  auto operator=(const reads_file_t &x) -> reads_file_t & = delete;
  auto operator=(reads_file_t &&) noexcept -> reads_file_t & = delete;
  reads_file_t(reads_file_t &&) noexcept = default;
  ~reads_file_t() noexcept = default;
  // clang-format on

  [[nodiscard]] friend auto
  is_active(const reads_file_t &x) -> bool {
    return x.self_->is_active_();
  }

  friend auto
  make_tasks(auto &reads_file, const std::int64_t n_chunks,
             const std::int32_t file_id, task_queue &tq,
             std::atomic_int32_t &n_tasks) -> void {
    reads_file.self_->make_tasks_(n_chunks, file_id, tq, n_tasks);
  }

private:
  // NOLINTNEXTLINE(cppcoreguidelines-special-member-functions)
  struct concept_t {
    // clang-format off
    virtual ~concept_t() = default;
    virtual auto make_tasks_(const std::int64_t, const std::int32_t,
                             task_queue &, std::atomic_int32_t &) -> void = 0;
    virtual auto is_active_() const -> bool = 0;
    // clang-format on
  };
  template <typename T> struct model : concept_t {
    explicit model(T x) : data_(std::move(x)) {}

    auto
    make_tasks_(const std::int64_t n_chunks, const std::int32_t file_id,
                task_queue &tq, std::atomic_int32_t &n_tasks) -> void override {
      make_tasks(data_, n_chunks, file_id, tq, n_tasks);
    }

    [[nodiscard]] auto
    is_active_() const -> bool override {
      return static_cast<bool>(data_);
    }

    T data_;
  };
  std::unique_ptr<concept_t> self_;
};

#endif  // SRC_READS_FILE_HPP_
