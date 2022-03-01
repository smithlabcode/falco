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

#include <string>
#include <cmath>

// Optional zlib usage
#include <zlib.h>

// Memory map
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/mman.h>

#ifdef USE_HTS
#include <htslib/sam.h>
#endif

#include "FalcoConfig.hpp"
#include "FastqStats.hpp"

/*************************************************************
 ******************** STREAM READER **************************
 *************************************************************/

// Generic class that does as much as possible without assuming file format
class StreamReader{
 public:
  // config on how to handle reads
  const bool do_sequence_hash,
             do_kmer,
             do_adapter,
             do_sliding_window,
             do_n_content,
             do_quality_base,
             do_sequence,
             do_gc_sequence,
             do_quality_sequence,
             do_tile,
             do_sequence_length;

  // This will tell me which character to look for to go to the next field
  const char field_separator;

  // This will tell me which character to look for to finish processing a record
  const char line_separator;

  // buffer size to store line 2 of each read statically
  const size_t buffer_size;

  /************ ADAPTER SEARCH ***********/
  const bool adapters_search_slow;
  const std::vector<std::string> adapter_seqs;

  const size_t num_adapters;
  const size_t adapter_size;
  const size_t adapter_mask;
  const std::array<size_t, Constants::max_adapters> adapters;

    // keep track of reads for which to do kmer and tile count
  static const size_t num_reads_for_tile = 10;
  static const size_t num_reads_for_kmer = 50;

  bool continue_storing_sequences;
  bool do_kmer_read;
  bool do_tile_read;

  size_t next_tile_read;
  size_t next_kmer_read;

  // Whether or not we have passed the buffer while reading and need to allocate
  // more space / use dynamically allocated space to process the base
  bool still_in_buffer;

  // Whether to just ignore per tile sequence quality if tiles are not in the
  // name
  bool tile_ignore;

  // Get a base from the sequence line
  char base_from_buffer;

  // the static buffer that takes the sequence data
  char *buffer;

  // pointer to current read character. Each file type should figure out how to
  // handle this as the file is processed
  char *cur_char;

  // Number of bases that have overflown the buffer
  size_t leftover_ind;

  /********* TILE PARSING ********/
  size_t num_colon;
  // tile value parsed from line 1 of each record
  size_t tile_cur;

  // the number of colons (:) needed to be seen until we know we are in a tile
  size_t tile_split_point;

  // Temp variables to be updated as you pass through the file
  size_t base_ind;  // 0,1,2 or 3
  size_t read_pos;  // which base we are at in the read
  size_t quality_value;  // to convert from ascii to number
  size_t cur_gc_count;  // Number of gc bases in read
  size_t cur_quality;  // Sum of quality values in read
  size_t num_bases_after_n;  // count of k-mers that reset at every N
  size_t cur_kmer;  // 32-mer hash as you pass through the sequence line
  size_t i; // general iterator

  // variables for gc model
  GCModelValue value;
  size_t truncated_length;
  size_t truncated_gc_count;
  size_t next_truncation;

  // Temporarily store line 2 out of 4 to know the base to which
  // quality characters are associated
  std::string leftover_buffer;
  std::string sequence_to_hash;  // sequence marked for duplication
  std::string filename;
  /************ FUNCTIONS TO PROCESS READS AND BASES ***********/
  // gets and puts bases from and to buffer
  inline void put_base_in_buffer();  // puts base in buffer or leftover
  inline void get_base_from_buffer();  // gets base from buffer or leftover

  // on the first tile-processed read, we will try to figure out how tiles
  // should be parsed
  inline void get_tile_split_position();
  inline void get_tile_value();

  inline void process_sequence_base_from_buffer(FastqStats &stats);
  inline void process_sequence_base_from_leftover(FastqStats &stats);
  inline void postprocess_sequence_line(FastqStats &stats);

  inline void process_quality_base_from_buffer(FastqStats &stats);
  inline void process_quality_base_from_leftover(FastqStats &stats);

  inline void postprocess_fastq_record(FastqStats &stats);

  /************ FUNCTIONS TO READ LINES IN DIFFERENT WAYS ***********/
  inline void read_fast_forward_line();  // run this to ignore a line
  inline void skip_separator();  // keep going forward while = separator
  inline void read_tile_line(FastqStats &stats);  // get tile from read name
  inline void read_sequence_line(FastqStats &stats);  // parse sequence
  inline void read_quality_line(FastqStats &stats);  // parse quality

  StreamReader(FalcoConfig &config, const size_t buffer_size,
               const char _field_separator, const char _line_separator);

  /************ FUNCTIONS TO IMPLEMENT BASED ON FILE FORMAT  ***********/
  virtual size_t load() = 0;
  virtual bool read_entry (FastqStats &stats, size_t &num_bytes_read) = 0;
  virtual bool is_eof() = 0;  // whether file has ended, for each file type
  virtual ~StreamReader() = 0;
};


/*******************************************************/
/*************** READ FASTQ RECORD *********************/
/*******************************************************/
class FastqReader : public StreamReader {
 private:
  static const size_t kChunkSize = (1<<20);
  char filebuf[kChunkSize];
  FILE *fileobj;

 public:
  FastqReader(FalcoConfig &fc, const size_t _buffer_size);

  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);
  ~FastqReader();
};

/*******************************************************/
/*************** READ FASTQ GZ RCORD *******************/
/*******************************************************/
class GzFastqReader : public StreamReader {
 private:
  static const size_t kChunkSize = (1<<20);
  char gzbuf[kChunkSize];
  gzFile fileobj;

 public:
  GzFastqReader(FalcoConfig &fc, const size_t _buffer_size);
  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);
  ~GzFastqReader();
};

/*******************************************************/
/*************** READ SAM RECORD ***********************/
/*******************************************************/

class SamReader : public StreamReader {
 private:
  // for uncompressed
  struct stat st;
  void *mmap_data;
  char *last;
  char *first;
 public:
  SamReader(FalcoConfig &fc, const size_t _buffer_size);
  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);
  ~SamReader();
};

#ifdef USE_HTS
/*******************************************************/
/*************** READ BAM RECORD ***********************/
/*******************************************************/

class BamReader : public StreamReader {
 private:
  htsFile *hts;
  bam_hdr_t *hdr;
  bam1_t *b;
  int rd_ret, fmt_ret;
  char *last;
 public:
  BamReader(FalcoConfig &fc, const size_t _buffer_size);
  size_t load();
  bool is_eof();
  bool read_entry(FastqStats &stats, size_t &num_bytes_read);
  ~BamReader();
};
#endif

#endif
