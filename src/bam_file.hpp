// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_BAM_FILE_HPP_
#define SRC_BAM_FILE_HPP_

#include "bamrec.hpp"
#include "falco_task.hpp"
#include "thread_pool_wrapper.hpp"  // for falco_thread_pool

#include <htslib/sam.h>

#include <algorithm>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <format>
#include <iterator>
#include <limits>
#include <memory>
#include <ranges>
#include <stdexcept>
#include <string>
#include <system_error>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

struct bam_buffer {
  static constexpr std::int32_t bam1_t_size = sizeof(bam1_t);
  static constexpr auto min_bytes_per_record = 256;
  std::int64_t n_recs{};
  std::int64_t n_bytes{};
  std::vector<bam1_t> recs;
  std::vector<std::uint8_t> data;
  // clang-format off
  explicit bam_buffer(const std::int64_t buf_size) :
    n_recs{buf_size / (bam1_t_size + min_bytes_per_record)},
    n_bytes{buf_size - n_recs * bam1_t_size},
    recs(n_recs), data(n_bytes, 0) {
    std::ranges::for_each(recs, [](auto &r) {
      bam_set_mempolicy(&r, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);
    });
  }
  [[nodiscard]] auto begin() {return std::begin(recs);}
  [[nodiscard]] auto begin() const {return std::cbegin(recs);}
  [[nodiscard]] auto end() {return std::end(recs); }
  [[nodiscard]] auto end() const {return std::cend(recs); }
  // clang-format on
  ~bam_buffer() {
    std::ranges::for_each(recs, [](auto &r) { bam_destroy1(&r); });
  }
};

struct bam_file {
  using rec_t = bamrec;
  static constexpr auto min_buf_size = static_cast<std::int64_t>(64 * 1024);
  static constexpr auto max_buf_size = std::numeric_limits<std::int64_t>::max();
  bam_buffer buf;
  std::unique_ptr<htsFile, int (*)(htsFile *)> f;
  std::unique_ptr<sam_hdr_t, void (*)(sam_hdr_t *)> h;
  bool hit_eof{};

  bam_file(const std::string &filename, const std::int64_t buf_size,
           falco_thread_pool &t) :
    buf(buf_size < min_buf_size ? min_buf_size : buf_size),
    f(hts_open(std::data(filename), "r"), &hts_close),
    h(sam_hdr_read(f.get()), &sam_hdr_destroy) {
    if (!f)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
    if (!h)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to read header: " + filename);
    if (t.n_threads() > 0) {
      // only use a thread pool if we have more than one thread
      const auto r = hts_set_thread_pool(f.get(), &t.t);
      if (r < 0)
        throw std::runtime_error("failed to set thread pool");
    }
  }

  // clang-format off
  // delete copy and assignment
  bam_file(const bam_file &) = delete;
  auto operator=(const bam_file &) -> bam_file & = delete;
  auto operator=(bam_file &&) noexcept -> bam_file & = delete;
  // default move for emplace
  bam_file(bam_file &&) noexcept = default;
  ~bam_file() noexcept = default;
  // clang-format on

  [[nodiscard]] operator bool() const { return !hit_eof; }

  auto
  load_next() -> const bam_file &;

  [[nodiscard]] auto
  get_chunks(std::int64_t n_chunks)
    -> std::vector<std::pair<rec_t::pos_t, rec_t::pos_t>>;
};

[[nodiscard]] auto
estimate_n_reads_bam(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t>;

[[nodiscard]] inline constexpr auto
is_active(const bam_file &reads_file) -> bool {
  return static_cast<bool>(reads_file);
}

[[nodiscard]] inline constexpr auto
inflate_only(bam_file &) -> bool {
  return false;
}

[[nodiscard]] inline constexpr auto
make_tasks_inflate(bam_file &) -> std::vector<task_t> {
  return {};
}

[[nodiscard]] auto
make_tasks(bam_file &reads_file,
           const std::int64_t n_chunks) -> std::vector<task_t>;

#endif  // SRC_BAM_FILE_HPP_
