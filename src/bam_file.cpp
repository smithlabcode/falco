// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bam_file.hpp"
#include "bamrec.hpp"
#include "falco_utils.hpp"

#include <htslib/hfile.h>
#include <htslib/sam.h>

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>

[[nodiscard]] auto
estimate_n_reads_bam(const std::string &filename)
  -> std::tuple<std::uint64_t, std::uint64_t, std::int64_t> {
  static constexpr auto max_n_reads = 128 * 1024;
  std::unique_ptr<htsFile, int (*)(htsFile *)> f(
    hts_open(std::data(filename), "r"), &hts_close);
  if (!f)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to open file: " + filename);
  std::unique_ptr<sam_hdr_t, void (*)(sam_hdr_t *)> h(sam_hdr_read(f.get()),
                                                      &sam_hdr_destroy);
  if (!h)
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "failed to read header: " + filename);

  const auto format = hts_get_format(f.get());
  if (!format)
    throw std::runtime_error("failed to identify file format: " + filename);

  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-union-access)
  const auto &fp = format->format == bam ? f->fp.bgzf->fp : f->fp.hfile;
  const auto pos_after_header = htell(fp);
  std::unique_ptr<bam1_t, void (*)(bam1_t *)> rec(bam_init1(), &bam_destroy1);
  std::uint64_t n_reads{};
  std::uint64_t total_read_len{};
  int r{};
  while (n_reads++ < max_n_reads &&
         (r = sam_read1(f.get(), h.get(), rec.get())) >= 0)
    total_read_len += rec->core.l_qseq;
  if (r < -1)  // error
    throw std::system_error(std::make_error_code(std::errc(errno)),
                            "error reading bam record from: " + filename);
  const auto pos_after_reads = htell(fp);
  const auto n_compressed_bytes = pos_after_reads - pos_after_header;
  const auto filesize = std::filesystem::file_size(filename);
  const auto estimate = static_cast<std::uint64_t>(
    as_frac(n_reads * (filesize - pos_after_header), n_compressed_bytes));
  return {estimate, total_read_len / n_reads, filesize};
}

auto
bam_file::load_next(const std::int32_t file_id, task_queue &tq,
                    std::atomic_int32_t &n_tasks) -> void {
  is_first_load = false;
  if (output_cursor > 0) {
    const auto n_to_keep = output_last - output_cursor;
    std::copy_n(std::cbegin(output_buffer) + output_cursor, n_to_keep,
                std::begin(input_buffer));
    input_last = n_to_keep;
  }
  br.reset();  // use like monotonic_buffer_resource
  auto in_itr = std::data(input_buffer) + input_last;
  while (br.task_ready() && has_in()) {
    auto task = br.get_decomp_task(in_itr);
    in_itr += task.size;
    ++n_tasks;
    tq.push(file_id, std::move(task));
  }
  input_last = std::distance(std::data(input_buffer), in_itr);
}

[[nodiscard]] inline auto
partition(auto itr,                     //
          const auto end,               //
          const std::int64_t n_chunks,  //
          const std::int32_t file_id,   //
          task_queue &tq,               //
          std::atomic_int32_t &n_tasks) {
  // ADS: this isn't working as desired: the end position of each part should be
  // the first record end past the 'end_itr' below unless end_itr == end
  const auto dist = std::distance(itr, end);
  const auto [chunk_size, remainder] = std::div(dist, n_chunks);
  bam_chunks_t positions;
  while (itr != end) {
    const auto dist = std::distance(itr, end);
    auto end_itr = itr + (dist < chunk_size ? dist : chunk_size);
    auto next_itr = bamrec::find_end_pos(itr, end_itr);
    if (next_itr == itr)
      break;
    ++n_tasks;
    tq.push(file_id, bam_task_t(itr, next_itr));
    itr = next_itr;
  }
  return itr;
}

auto
bam_file::get_chunks(const std::int64_t n_chunks,  //
                     const std::int32_t file_id,   //
                     task_queue &tq,               //
                     std::atomic_int32_t &n_tasks) -> void {
  std::swap(input_buffer, output_buffer);
  std::swap(output_last, input_last);
  auto itr = std::data(output_buffer);
  const auto end = itr + output_last;
  if (bh)  // if we are still parsing the header
    itr = bh.update(itr, end);
  if (std::distance(itr, end) > 0)
    itr = partition(itr, end, n_chunks, file_id, tq, n_tasks);
  output_cursor = std::distance(std::data(output_buffer), itr);
}

auto
bam_file::make_tasks(const std::int64_t n_chunks, const std::int32_t file_id,
                     task_queue &tq, std::atomic_int32_t &n_tasks) -> void {
  get_chunks(n_chunks, file_id, tq, n_tasks);
  had_last_chunks = !static_cast<bool>(br);  // wait until now to set this?
}
