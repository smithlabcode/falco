/* Copyright (C) 2019 Guilherme De Sena Brandine and
 *                    Andrew D. Smith
 * Authors: Guilherme De Sena Brandine, Andrew Smith
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

#ifndef STREAMREADER_HPP
#define STREAMREADER_HPP

#include <bitset>
#include <cmath>
#include <string>

#include <sys/stat.h>
#include <zlib.h>  // Optional zlib usage

#ifdef USE_HTS
#ifdef USE_HTS
#include "bamxx/bamxx.hpp"
#endif
// #include <htslib/bgzf.h>
// #include <htslib/sam.h>
#endif

#include "FalcoConfig.hpp"
#include "FastqStats.hpp"

/*************************************************************
 ******************** STREAM READER **************************
 *************************************************************/

// Generic struct that does as much as possible without assuming file format
struct StreamReader {
  // keep track of reads for which to do kmer and tile count
  static constexpr std::size_t num_reads_for_tile = 10;
  static constexpr std::size_t num_reads_for_kmer = 50;
  static constexpr std::size_t check_bytes_read_mask = 65535;

  // config on how to handle reads
  const bool do_sequence_hash{};
  const bool do_kmer{};
  const bool do_adapter{};
  const bool do_adapter_optimized{};
  const bool do_sliding_window{};
  const bool do_n_content{};
  const bool do_quality_base{};
  const bool do_sequence{};
  const bool do_gc_sequence{};
  const bool do_quality_sequence{};
  const bool do_tile{};
  const bool do_sequence_length{};

  // This will tell me which character to look for to go to the next field
  const char field_separator{};

  // This will tell me which character to look for to finish processing a record
  const char line_separator{};

  // buffer size to store line 2 of each read statically
  const std::size_t buffer_size{};
  const std::size_t read_step{};
  // the number of colons (:) needed to be seen until we know we are in a tile
  const std::size_t tile_split_point{};

  // Whether to just ignore per tile sequence quality if tiles are not in the
  // name
  const bool tile_ignore{};

  /************ ADAPTER SEARCH ***********/
  const bool do_adapters_slow{};
  const std::vector<std::string> adapter_seqs;

  const std::size_t num_adapters{};
  const std::size_t adapter_size{};
  const std::size_t adapter_mask{};
  const std::array<std::size_t, Constants::max_adapters> adapters{};
  std::bitset<Constants::max_adapters> adapters_found{};

  const std::string filename;

  bool continue_storing_sequences{};
  bool do_read{};
  bool do_kmer_read{};
  bool do_tile_read{};

  std::size_t next_read{};
  std::size_t next_tile_read{};
  std::size_t next_kmer_read{};

  // Whether or not we have passed the buffer while reading and need to allocate
  // more space / use dynamically allocated space to process the base
  bool still_in_buffer{};

  // Get a base from the sequence line
  char base_from_buffer{};

  // the static buffer that takes the sequence data
  char *buffer{};

  // pointer to current read character. Each file type should figure out how to
  // handle this as the file is processed
  char *cur_char{};

  // Number of bases that have overflown the buffer
  std::size_t leftover_ind{};

  /********* TILE PARSING ********/
  // tile value parsed from line 1 of each record
  std::size_t tile_cur{};

  // Temp variables to be updated as you pass through the file
  std::size_t read_pos{};           // which base we are at in the read
  std::size_t quality_value{};      // to convert from ascii to number
  std::size_t cur_gc_count{};       // Number of gc bases in read
  std::size_t cur_quality{};        // Sum of quality values in read
  std::size_t num_bases_after_n{};  // count of k-mers that reset at every N
  std::size_t cur_kmer{};  // 32-mer hash as you pass through the sequence line

  // variables for gc model
  GCModelValue value;
  std::size_t truncated_length{};
  std::size_t truncated_gc_count{};
  std::size_t next_truncation{};

  // Temporarily store line 2 out of 4 to know the base to which quality
  // characters are associated
  std::string leftover_buffer;
  std::string sequence_to_hash;  // sequence marked for duplication

  StreamReader(FalcoConfig &config, const std::size_t buffer_size,
               const char field_separator, const char line_separator);

  /************ FUNCTIONS TO PROCESS READS AND BASES ***********/
  // gets and puts bases from and to buffer
  void
  put_base_in_buffer();  // puts base in buffer or leftover

  void
  get_base_from_buffer();  // gets base from buffer or leftover

  // on the first tile-processed read, we will try to figure out how tiles
  // should be parsed
  void
  get_tile_value();

  void
  process_sequence_base_from_buffer(FastqStats &stats);

  void
  process_sequence_base_from_leftover(FastqStats &stats);

  void
  postprocess_sequence_line(FastqStats &stats);

  void
  process_quality_base_from_buffer(FastqStats &stats);

  void
  process_quality_base_from_leftover(FastqStats &stats);

  void
  postprocess_fastq_record(FastqStats &stats);

  /************ FUNCTIONS TO READ LINES IN DIFFERENT WAYS ***********/
  void
  read_fast_forward_line();  // run this to ignore a line

  void
  read_fast_forward_line_eof();  // run this to ignore a line until EOF

  void
  skip_separator();  // keep going forward while = separator

  void
  read_tile_line(FastqStats &stats);  // get tile from read name

  void
  read_sequence_line(FastqStats &stats);  // parse sequence

  void
  read_quality_line(FastqStats &stats);  // parse quality

  /************ FUNCTIONS FOR PROGRESS BAR ***********/
  bool
  check_bytes_read(const std::size_t line_num);

  /************ FUNCTIONS TO IMPLEMENT BASED ON FILE FORMAT  ***********/
  virtual std::size_t
  load() = 0;

  virtual bool
  read_entry(FastqStats &stats, std::size_t &num_bytes_read) = 0;

  virtual bool
  is_eof() = 0;  // whether file has ended, for each file type

  virtual ~StreamReader() = 0;
};

struct FastqReader : public StreamReader {
  static constexpr std::size_t RESERVE_SIZE = (1 << 26);
  std::vector<char> file_buffer;
  std::FILE *fileobj{};

  FastqReader(FalcoConfig &fc, const std::size_t _buffer_size);

  std::size_t
  load();

  bool
  is_eof();

  bool
  read_entry(FastqStats &stats, std::size_t &num_bytes_read);

  ~FastqReader();
};

struct GzFastqReader : public StreamReader {
  static constexpr std::size_t RESERVE_SIZE = (1 << 26);
  std::vector<char> gzbuf;
  gzFile fileobj;

  GzFastqReader(FalcoConfig &fc, const std::size_t _buffer_size);

  std::size_t
  load();

  bool
  is_eof();

  bool
  read_entry(FastqStats &stats, std::size_t &num_bytes_read);

  ~GzFastqReader();
};

struct SamReader : public StreamReader {
  static constexpr std::size_t RESERVE_SIZE = (1 << 26);
  std::vector<char> file_buffer;
  std::FILE *fileobj{};

  SamReader(FalcoConfig &fc, const std::size_t _buffer_size);

  std::size_t
  load();

  bool
  is_eof();

  bool
  read_entry(FastqStats &stats, std::size_t &num_bytes_read);

  ~SamReader();
};

#endif
