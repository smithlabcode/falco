// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_READS_FILE_HPP_
#define SRC_READS_FILE_HPP_

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <variant>
#include <vector>

class reads_file_t {
public:
  // clang-format off
  template <typename T> reads_file_t(T x) : self_(new model<T>(std::move(x))) {}
  reads_file_t(const reads_file_t &x) = delete;
  auto operator=(const reads_file_t &x) -> reads_file_t & = delete;
  auto operator=(reads_file_t &&) noexcept -> reads_file_t & = delete;
  reads_file_t(reads_file_t &&) noexcept = default;
  // clang-format on

  [[nodiscard]] friend auto
  is_active(const reads_file_t &x) -> bool {
    return x.self_->is_active_();
  }

  [[nodiscard]] friend auto
  make_tasks(auto &reads_file,
             const std::int64_t n_chunks) -> std::vector<task_t> {
    return reads_file.self_->make_tasks_(n_chunks);
  }

  [[nodiscard]] friend auto
  make_tasks_inflate(auto &reads_file) -> std::vector<task_t> {
    return reads_file.self_->make_tasks_inflate_();
  }

  [[nodiscard]] friend auto
  inflate_only(reads_file_t &x) -> bool {
    return x.self_->inflate_only_();
  }

  friend auto
  read_data(reads_file_t &x) -> void {
    x.self_->read_data_();
  }

private:
  struct concept_t {
    // clang-format off
    virtual ~concept_t() = default;
    virtual auto make_tasks_(const std::int64_t) -> std::vector<task_t> = 0;
    virtual auto make_tasks_inflate_() -> std::vector<task_t> = 0;
    virtual auto inflate_only_() -> bool = 0;
    virtual auto is_active_() const -> bool = 0;
    virtual auto read_data_() -> void = 0;
    // clang-format on
  };
  template <typename T> struct model : concept_t {
    model(T x) : data_(std::move(x)) {}

    [[nodiscard]] auto
    make_tasks_(const std::int64_t n_chunks) -> std::vector<task_t> override {
      return make_tasks(data_, n_chunks);
    }

    [[nodiscard]] auto
    make_tasks_inflate_() -> std::vector<task_t> override {
      return make_tasks_inflate(data_);
    }

    [[nodiscard]] auto
    is_active_() const -> bool override {
      return static_cast<bool>(data_);
    }

    auto
    inflate_only_() -> bool override {
      return inflate_only(data_);
    }

    auto
    read_data_() -> void override {
      read_data(data_);
    }

    T data_;
  };
  std::unique_ptr<concept_t> self_;
};

#endif  // SRC_READS_FILE_HPP_
