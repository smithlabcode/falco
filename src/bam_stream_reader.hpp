/* Copyright (C) 2026 Andrew D. Smith
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 */

#ifndef BAM_STREAM_READER_HPP_
#define BAM_STREAM_READER_HPP_

#include "StreamReader.hpp"

#include "bamxx/bamxx.hpp"

struct BamReader : public StreamReader {
  bamxx::bam_tpool tp;
  bamxx::bam_in hts;
  bamxx::bam_header hdr;
  bamxx::bam_rec aln;
  int rd_ret{};
  int fmt_ret{};
  char *last{};

  BamReader(FalcoConfig &fc, const std::size_t _buffer_size);

  [[nodiscard]] std::size_t
  load();

  [[nodiscard]] bool
  is_eof();

  [[nodiscard]] bool
  read_entry(FastqStats &stats, std::size_t &num_bytes_read);

  // Specially made for BamReader to work directly with bam1_t
  void
  read_sequence_line(FastqStats &stats);

  void
  read_quality_line(FastqStats &stats);  // parse quality

  void
  put_base_in_buffer(const std::size_t pos);  // puts base in buffer or leftover

  ~BamReader() {}

  void
  process_sequence_base_from_buffer(FastqStats &stats);

  void
  postprocess_fastq_record(FastqStats &stats);
};

#endif  // BAM_STREAM_READER_HPP_
