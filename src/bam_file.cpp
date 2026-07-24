// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "bam_file.hpp"
#include "falco_utils.hpp"

#include <htslib/bgzf.h>   // for BGZF
#include <htslib/hfile.h>  // for htell

#include <cstdint>
#include <filesystem>
#include <memory>
#include <string>
#include <system_error>

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
bam_file::load_next() -> const bam_file & {
  // ADS: need to make sure the buffer starts at the proper alignment
  const auto align = [](const auto l) {
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-magic-numbers)
    return static_cast<std::uint32_t>(l + 7) & ~7u;
  };
  auto &recs = buf.recs;
  // release the final BAM record in the previous load_next, if any
  const auto prev_n = buf.n_recs;
  if (prev_n > 1 &&
      (bam_get_mempolicy(&recs[prev_n - 1]) & BAM_USER_OWNS_DATA) == 0)
    bam_destroy1(&recs[prev_n - 1]);

  auto n_bytes = 0u;
  auto n_recs = 0u;

  // ADS: if data buffer capacity exceeded, BAM_USER_OWNS_DATA check fails and
  // loop terminates; one record will allocate space for data from the heap
  while (n_recs < std::size(recs)) {
    auto &rec = recs[n_recs];
    bam_set_mempolicy(&rec, BAM_USER_OWNS_STRUCT | BAM_USER_OWNS_DATA);
    // NOLINTNEXTLINE(cppcoreguidelines-pro-bounds-pointer-arithmetic)
    rec.data = std::data(buf.data) + n_bytes;
    rec.m_data = std::size(buf.data) - n_bytes;
    const auto r = sam_read1(f.get(), h.get(), &rec);
    if (r < -1)
      throw std::runtime_error("error reading bam file");
    if (r == -1) {
      hit_eof = true;
      break;
    }
    ++n_recs;
    // if there is no space for the record, stop
    if ((bam_get_mempolicy(&rec) & BAM_USER_OWNS_DATA) == 0)
      break;  // no more space
    // round up to 8 bytes for memory alignment
    rec.m_data = align(rec.l_data);
    n_bytes += rec.m_data;
  }
  buf.n_recs = n_recs;
  buf.n_bytes = n_bytes;
  return *this;
}

[[nodiscard]] auto
bam_file::get_chunks(std::int64_t n_chunks) -> bam_chunks_t {
  if (buf.n_recs == 0)
    return bam_chunks_t(1, {std::cbegin(buf), std::cbegin(buf)});
  n_chunks = std::min(n_chunks, buf.n_recs);
  const auto [chunk_size, remainder] = std::div(buf.n_recs, n_chunks);
  const auto buffer = std::cbegin(buf);
  std::int64_t start_pos{};
  bam_chunks_t chunks(n_chunks);
  for (const auto chunk_idx : std::views::iota(0, n_chunks)) {
    const auto stop_pos = start_pos + chunk_size + (chunk_idx < remainder);
    chunks[chunk_idx] = {buffer + start_pos, buffer + stop_pos};
    start_pos = stop_pos;
  }
  return chunks;
}

[[nodiscard]] auto
make_tasks(bam_file &reads_file,
           const std::int64_t n_chunks) -> std::vector<task_t> {
  reads_file.load_next();
  auto chunks = reads_file.get_chunks(n_chunks);
  std::vector<task_t> tasks;
  tasks.reserve(std::size(chunks));
  for (const auto &chunk : chunks)
    tasks.push_back(bam_task_t(chunk.first, chunk.second));
  return tasks;
}
