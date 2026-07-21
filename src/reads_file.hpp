// SPDX-License-Identifier: MIT; (c) 2026 Andrew D Smith

#ifndef SRC_READS_FILE_HPP_
#define SRC_READS_FILE_HPP_

#include "bam_file.hpp"
#include "bgzf_block.hpp"
#include "fastq_file.hpp"
#include "fqrec.hpp"

#include <cstdint>
#include <cstdlib>
#include <memory>
#include <variant>
#include <vector>

struct fq_task_t {
  fqrec::pos_t beg{};
  fqrec::pos_t end{};
};

struct bam_task_t {
  bamrec::pos_t beg{};
  bamrec::pos_t end{};
};

using task_t =
  std::variant<std::monostate, fq_task_t, bam_task_t, bgzf_block_t>;

[[nodiscard]] inline auto
make_tasks_impl(bam_chunks_t chunks) -> std::vector<task_t> {
  std::vector<task_t> tasks;
  tasks.reserve(std::size(chunks));
  for (const auto &chunk : chunks)
    tasks.push_back(bam_task_t(chunk.first, chunk.second));
  return tasks;
}

[[nodiscard]] inline auto
make_tasks_impl(fq_chunks_t chunks) -> std::vector<task_t> {
  std::vector<task_t> tasks;
  tasks.reserve(std::size(chunks));
  for (const auto &chunk : chunks)
    tasks.push_back(fq_task_t(chunk.first, chunk.second));
  return tasks;
}

using bgzf_chunks_t = std::vector<bgzf_block_t>;

[[nodiscard]] inline auto
inflate_only(fastq_bgzf_file &f) -> bool {
  return f.is_first_load;
}

[[nodiscard]] inline auto
make_tasks_inflate_impl(bgzf_chunks_t chunks) -> std::vector<task_t> {
  std::vector<task_t> tasks;
  tasks.reserve(std::size(chunks));
  for (auto &chunk : chunks)
    tasks.emplace_back(std::move(chunk));
  return tasks;
}

[[nodiscard]] inline auto
make_tasks(fastq_bgzf_file &reads_file,
           const std::int64_t n_chunks) -> std::vector<task_t> {
  const auto chunks = reads_file.get_chunks(n_chunks);
  reads_file.had_last_chunks = reads_file.hit_eof;
  return make_tasks_impl(std::move(chunks));
}

[[nodiscard]] inline auto
make_tasks_inflate(fastq_bgzf_file &reads_file) -> std::vector<task_t> {
  reads_file.is_first_load = false;
  return make_tasks_inflate_impl(reads_file.load_next());
}

template <typename T>
concept concrete_reads_file_t =        //
  std::same_as<T, fastq_file> ||       //
  std::same_as<T, fastq_bgzf_file> ||  //
  std::same_as<T, fastq_gz_file> ||    //
  std::same_as<T, bam_file>;

[[nodiscard]] inline auto
inflate_only([[maybe_unused]] concrete_reads_file_t auto &f) -> bool {
  return false;
}

[[nodiscard]] inline auto
make_tasks_inflate([[maybe_unused]] concrete_reads_file_t auto &reads_file)
  -> std::vector<task_t> {
  return {};
}

[[nodiscard]] inline auto
make_tasks(concrete_reads_file_t auto &reads_file,
           const std::int64_t n_chunks) -> std::vector<task_t> {
  reads_file.load_next();
  return make_tasks_impl(reads_file.get_chunks(n_chunks));
}

[[nodiscard]] inline auto
is_active(const concrete_reads_file_t auto &reads_file) -> bool {
  return static_cast<bool>(reads_file);
}

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

  friend auto
  inflate_only(reads_file_t &x) -> bool {
    return x.self_->inflate_only_();
  }

private:
  struct concept_t {
    // clang-format off
    virtual ~concept_t() = default;
    virtual auto make_tasks_(const std::int64_t) -> std::vector<task_t> = 0;
    virtual auto make_tasks_inflate_() -> std::vector<task_t> = 0;
    virtual auto inflate_only_() -> bool = 0;
    virtual auto is_active_() const -> bool = 0;
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
      return is_active(data_);
    }

    auto
    inflate_only_() -> bool override {
      return inflate_only(data_);
    }

    T data_;
  };
  std::unique_ptr<concept_t> self_;
};

#endif  // SRC_READS_FILE_HPP_
