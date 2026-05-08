/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SRC_BAM_FILE_HPP_
#define SRC_BAM_FILE_HPP_

#include "falco_utils.hpp"

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <cstdint>
#include <cstring>
#include <filesystem>
#include <limits>
#include <memory>
#include <new>
#include <numeric>
#include <ranges>
#include <string>
#include <system_error>
#include <thread>
#include <utility>
#include <vector>

template <class T> struct aligned_allocator {
  static constexpr std::size_t align_at = 8;
  typedef T value_type;
  [[nodiscard]] auto
  allocate(const std::size_t n) -> T * {
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(T))
      throw std::bad_array_new_length();
    if (auto p = static_cast<T *>(std::aligned_alloc(align_at, n * sizeof(T))))
      return p;
    throw std::bad_alloc();
  }
  auto
  deallocate(T *p, [[maybe_unused]] const std::size_t n) noexcept {
    std::free(p);
  }
};

// clang-format off
struct bam_seq_itr {
  // ADS: this class exists to model a char * as is used for the read sequence
  // in the case of FASTQ files. It doesn't have all the functions usually
  // defined for iterators, but this is by design.
  const std::uint8_t *p{};
  std::uint64_t i{};
  explicit bam_seq_itr(const auto p) : p{p} {}
  [[nodiscard]] auto operator<=>(const bam_seq_itr &) const = default;
  [[nodiscard]] auto operator*() const { return seq_nt16_str[bam_seqi(p, i)]; }
  [[nodiscard]] auto operator[](const auto j) const { return seq_nt16_str[bam_seqi(p, i + j)]; }
  auto operator+(const auto rhs) {
    bam_seq_itr tmp(*this);
    tmp.i += rhs;
    return tmp;
  }
  auto operator++(int) {
    bam_seq_itr tmp(*this);
    ++i;
    return tmp;
  }
  auto operator++() {
    ++i;
    return *this;
  }
};

struct bamrec {
  using pos_t = std::vector<bam1_t>::const_iterator;
  std::int32_t l_qname{};
  std::int32_t l_qseq{};
  const char *n{};
  const std::uint8_t *r{};
  const std::uint8_t *q{};
  bamrec() = default;
  explicit bamrec(const bam1_t &b) :
    l_qname{b.core.l_qname},
    l_qseq{b.core.l_qseq},
    n{bam_get_qname(&b)},
    r{bam_get_seq(&b)},
    q{bam_get_qual(&b)}
  {}
  [[nodiscard]] operator bool() const { return n != nullptr; }
};

[[nodiscard]] constexpr auto get_name(const bamrec &rec) { return rec.n; }
[[nodiscard]] constexpr auto get_name_end(const bamrec &rec) { return rec.n + rec.l_qname; }

[[nodiscard]] constexpr auto get_seq(const bamrec &rec) { return bam_seq_itr(rec.r); }
[[nodiscard]] constexpr auto get_seq_end(const bamrec &rec) { return bam_seq_itr(rec.r) + rec.l_qseq; }
[[nodiscard]] constexpr auto get_seq_size(const bamrec &rec) { return rec.l_qseq; }

[[nodiscard]] constexpr auto get_qual(const bamrec &rec) { return rec.q; }
[[nodiscard]] constexpr auto get_qual_end(const bamrec &rec) { return *rec.q == 0xff ? rec.q : rec.q + rec.l_qseq; }
[[nodiscard]] constexpr auto get_qual_size(const bamrec &rec) { return rec.l_qseq; }
// clang-format on

[[nodiscard]] auto
to_string(const bamrec &b) {
  // converts the bam record to FASTQ format for visualization
  static constexpr auto quality_score_offset = 33;
  static constexpr auto other = "\n+\n";
  static constexpr auto qual_missing_code = 0xff;  // from sam.c
  auto name = std::format("@{}\n", std::string(b.n, b.l_qname));
  std::string read;
  auto seq_itr = get_seq(b);
  while (seq_itr != get_seq_end(b))
    read += *seq_itr++;
  std::string qual;
  auto qual_itr = get_qual(b);
  if (*qual_itr == qual_missing_code)
    qual = std::string(get_qual_size(b), 'B');
  else
    while (qual_itr != get_qual_end(b))
      qual += (quality_score_offset + *qual_itr++);
  return name + read + other + qual + '\n';
}

static inline constexpr std::int32_t bam1_t_size = sizeof(bam1_t);

struct bam_buffer {
  static constexpr auto min_bytes_per_record = 256;
  std::int32_t n_recs{};
  std::int32_t n_bytes{};
  std::vector<bam1_t> recs;
  // ADS: ensure alignment within the vector; I can't produce an allocation that
  // is not aligned to 8 bytes, but I can't find a reason to assume this happens
  std::vector<std::uint8_t, aligned_allocator<std::uint8_t>> data;
  // clang-format off
  explicit bam_buffer(const std::int32_t buf_size) :
    n_recs{buf_size / (bam1_t_size + min_bytes_per_record)},
    n_bytes{buf_size - n_recs * bam1_t_size},
    recs(n_recs),
    data(n_bytes, 0)
  {}
  [[nodiscard]] auto begin() {return std::begin(recs);}
  [[nodiscard]] auto begin() const {return std::cbegin(recs);}
  [[nodiscard]] auto end() {return std::end(recs); }
  [[nodiscard]] auto end() const {return std::cend(recs); }
  // clang-format on
};

[[nodiscard]] inline auto
get_next(bamrec::pos_t &cursor,
         [[maybe_unused]] const bamrec::pos_t end_itr) -> bamrec {
  auto tmp = cursor;
  ++cursor;
  return bamrec(*tmp);
}

struct bam_file {
  using rec_t = bamrec;
  static constexpr auto min_buf_size = 16 * 4096;
  static constexpr auto max_buf_size = std::numeric_limits<std::int32_t>::max();
  bam_buffer buf;
  falco_thread_pool t;
  std::unique_ptr<htsFile, int (*)(htsFile *)> f;
  std::unique_ptr<sam_hdr_t, void (*)(sam_hdr_t *)> h;
  bool hit_eof{};

  bam_file(const std::string &filename, const std::int32_t buf_size,
           const std::uint32_t n_threads = 1) :
    buf(std::clamp(buf_size, min_buf_size, max_buf_size)), t(n_threads),
    f(hts_open(std::data(filename), "r"), &hts_close),
    h(sam_hdr_read(f.get()), &sam_hdr_destroy) {
    if (!f)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to open file: " + filename);
    if (!h)
      throw std::system_error(std::make_error_code(std::errc(errno)),
                              "failed to read header: " + filename);
    if (n_threads > 1) {
      // only use a thread pool if we have more than one thread
      const auto r = hts_set_thread_pool(f.get(), &t.t);
      if (r < 0)
        throw std::runtime_error("failed to set thread pool");
    }
  }

  [[nodiscard]] operator bool() const { return !hit_eof; }

  auto
  load_next() -> const bam_file & {
    // ADS: need to make sure the buffer starts at the proper alignment
    const auto align = [](const auto l) {
      return static_cast<std::uint32_t>(l + 7) & ~7u;
    };
    auto n_bytes = 0u;
    auto n_recs = 0u;
    // ADS: if data buffer capacity exceeded, BAM_USER_OWNS_DATA check fails and
    // loop terminates
    while (n_recs < std::size(buf.recs)) {
      auto &rec = buf.recs[n_recs];
      bam_set_mempolicy(&rec, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);
      rec.data = std::data(buf.data) + n_bytes;
      rec.m_data = std::size(buf.data) - n_bytes;
      const auto r = sam_read1(f.get(), h.get(), &rec);
      if (r < -1)
        throw std::runtime_error("error reading bam file");
      if (r == -1) {
        hit_eof = true;
        break;
      }
      // if there is no space for the record, stop
      if ((bam_get_mempolicy(&rec) & BAM_USER_OWNS_DATA) == 0) {
        ++n_recs;  // include last record in count
        break;     // no more space
      }
      // round up to 8 bytes for memory alignment
      rec.m_data = align(rec.l_data);
      n_bytes += rec.m_data;
      ++n_recs;
    }
    buf.n_recs = n_recs;
    buf.n_bytes = n_bytes;
    return *this;
  }
};

[[nodiscard]] static inline auto
get_chunks(const bam_file &bf, const std::int32_t n_chunks)
  -> std::vector<std::pair<bamrec::pos_t, bamrec::pos_t>> {
  const auto [chunk_size, remainder] = std::div(bf.buf.n_recs, n_chunks);
  const auto buffer = std::cbegin(bf.buf);
  auto start_pos = 0;
  std::vector<std::pair<bamrec::pos_t, bamrec::pos_t>> chunks(n_chunks);
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    chunks[chunk_idx] = {buffer + start_pos, buffer + stop_pos};
    start_pos = stop_pos;
  }
  return chunks;
}

#endif  // SRC_BAM_FILE_HPP_
