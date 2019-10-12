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

#include "StreamReader.hpp"

#include <vector>
#include <cstring>

using std::string;
using std::vector;
using std::runtime_error;

/****************************************************/
/***************** STREAMREADER *********************/
/****************************************************/
StreamReader::StreamReader(FalcoConfig &config,
                           const size_t _buffer_size,
                           const char _field_separator,
                           const char _line_separator) :
  // I have to pass the config skips as const to read them fast
  do_duplication(config.do_duplication),
  do_kmer(config.do_kmer),
  do_n_content(config.do_n_content),
  do_overrepresented(config.do_overrepresented),
  do_quality_base(config.do_quality_base),
  do_sequence(config.do_sequence),
  do_gc_sequence(config.do_gc_sequence),
  do_quality_sequence(config.do_quality_sequence),
  do_tile(config.do_tile),
  do_adapter(config.do_adapter),
  do_sequence_length(config.do_sequence_length),
  do_nogroup(config.nogroup),

  // Here are the usual stream reader configs
  buffer_size(_buffer_size),
  field_separator(_field_separator),
  line_separator(_line_separator) {

  // Allocates buffer to temporarily store reads
  buffer = new char[buffer_size + 1];
  buffer[buffer_size] = '\0';

  // duplication init
  continue_storing_sequences = true;

  // Tile init
  tile_ignore = !do_tile;  // early ignore tile if asked to skip it
  tile_cur = 0;
  tile_split_point = 0;
  next_tile_read = 0;

  // kmer init
  next_kmer_read = 0;

  // Subclasses will use this to deflate if necessary
  filename = config.filename;
}

// Makes sure that any subclass deletes the buffer
StreamReader::~StreamReader() {
  delete buffer;
}

/*******************************************************/
/*************** BUFFER MANAGEMENT *********************/
/*******************************************************/
// puts base either on buffer or leftover
inline void
StreamReader::put_base_in_buffer() {
  base_from_buffer = *cur_char;
  if (write_to_buffer) {
    buffer[read_pos] = base_from_buffer;
  }
  else {
    if (leftover_ind == leftover_buffer.size())
      leftover_buffer.push_back(base_from_buffer);
    else
      leftover_buffer[leftover_ind] = base_from_buffer;
  }
}

// Gets base from either buffer or leftover
inline void
StreamReader::get_base_from_buffer() {
  if (read_from_buffer) {
    base_from_buffer = buffer[read_pos];
  }
  else {
    base_from_buffer = leftover_buffer[leftover_ind];
  }
}

/*******************************************************/
/*************** FAST FOWARD ***************************/
/*******************************************************/

// Keeps going forward while the current character is a separator
inline void
StreamReader::skip_separator() {
  for (; *cur_char == field_separator; ++cur_char) {}
}

// Skips lines that are not relevant
inline void
StreamReader::read_fast_forward_line() {
  for (; *cur_char != field_separator; ++cur_char) {}
}

/*******************************************************/
/*************** TILE PROCESSING ***********************/
/*******************************************************/
// Parse the comment
inline void
StreamReader::get_tile_split_position() {
  num_colon = 0;

  // Count colons to know the formatting pattern
  for (; *cur_char != field_separator; ++cur_char) {
    if (*cur_char == ':') ++num_colon;
  }

  // Copied from fastqc
  if (num_colon >= 6) {
    tile_split_point = 4;
  }
  else if (num_colon >=4) {
    tile_split_point = 2;
  }

  // We will not get a tile out of this
  else {
    tile_ignore = true;
  }
}

inline void
StreamReader::get_tile_value() {
  tile_cur = 0;
  num_colon = 0;
  for (; *cur_char != field_separator; ++cur_char) {
    if (*cur_char == ':') ++num_colon;
    if (num_colon == tile_split_point) {
      ++cur_char;

      // parse till next colon or \n
      for (; (*cur_char != ':') && (*cur_char != field_separator); ++cur_char) {
        tile_cur = tile_cur*10 + (*cur_char - '0');
      }
      ++num_colon;
    }
  }
}

// Gets the tile from the sequence name (if applicable)
inline void
StreamReader::read_tile_line(FastqStats &stats) {
  // if there is no tile information in the fastq header, fast
  // forward this line
  if (tile_ignore) {
    read_fast_forward_line();
    return;
  }

  // Fast forward if this is not a tile line
  if (stats.num_reads != next_tile_read) {
    read_fast_forward_line();
    return;
  }

  // We haven't parsed the first line to know the split point
  if (tile_split_point == 0)  {
    get_tile_split_position();
  }
  else {
    get_tile_value();
    // allocate vector for tile if it doesn't exist
    if (stats.tile_position_quality.count(tile_cur) == 0) {
      stats.tile_count[tile_cur] = 0;
      stats.tile_position_quality[tile_cur] =
        vector<double> (stats.max_read_length, 0.0);
      stats.tile_position_count[tile_cur] =
        vector<size_t> (stats.max_read_length, 0);
    }
  }
}


/*******************************************************/
/*************** SEQUENCE PROCESSING *******************/
/*******************************************************/

// This is probably the most important function for speed, so it must be really
// optimized at all times
inline void
StreamReader::process_sequence_base_from_buffer(FastqStats &stats) {
  // I will count the Ns even if asked to ignore, as checking ifs take time
  if (base_from_buffer == 'N') {
    stats.n_base_count[read_pos]++;
    num_bases_after_n = 1;  // start over the current kmer
  }

  // ATGC bases
  else {
    // two bit base index
    base_ind = actg_to_2bit(base_from_buffer);

    // increments basic statistic counts
    cur_gc_count += (base_ind & 1);
    stats.base_count[(read_pos << stats.kBitShiftNucleotide) | base_ind]++;

    if (do_adapter) {
      if (stats.num_reads == next_kmer_read) {
        if (read_pos < stats.kKmerMaxBases) {
          cur_kmer = ((cur_kmer << stats.kBitShiftNucleotide) | base_ind);

          // registers k-mer if seen at least k nucleotides since the last n
          if (num_bases_after_n == stats.kmer_size) {
            stats.kmer_count[(read_pos << stats.kBitShiftKmer)
                             | (cur_kmer & stats.kmer_mask)]++;
          }

          else {
            num_bases_after_n++;
          }
        }
      }
    }
  }
}

// slower version of process_sequence_base_from_buffer that dynamically
// allocates
inline void
StreamReader::process_sequence_base_from_leftover(FastqStats &stats) {
  if (base_from_buffer == 'N') {
    stats.long_n_base_count[leftover_ind]++;
    num_bases_after_n = 1;  // start over the current kmer
  }

  // ATGC bases
  else {
    // two bit base index
    base_ind = actg_to_2bit(base_from_buffer);

    // increments basic statistic counts
    cur_gc_count += (base_ind & 1);
    stats.long_base_count[(leftover_ind << stats.kBitShiftNucleotide)
                          | base_ind]++;

    // WE WILL NOT DO KMER STATS OUTSIDE OF BUFFER
  }
}

// Gets statistics after reading the entire sequence line
inline void
StreamReader::postprocess_sequence_line(FastqStats &stats) {
  if (do_sequence_length) {
    // read length frequency histogram
    if ((read_pos != 0) && (read_pos <= stats.kNumBases))
      stats.read_length_freq[read_pos - 1]++;
    else
      stats.long_read_length_freq[leftover_ind - 1]++;
  }

  // Updates maximum read length if applicable
  if (read_pos > stats.max_read_length) {
    stats.max_read_length = read_pos;
  }

  // Registers GC % in the bin truncated to the nearest integer
  if (do_gc_sequence) {
    stats.total_gc += cur_gc_count;
    stats.gc_count[round(100 * cur_gc_count / static_cast<double>(read_pos))]++;
  }
}

// Reads the line that has the biological sequence
inline void
StreamReader::read_sequence_line(FastqStats &stats) {
  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  num_bases_after_n = 1;
  write_to_buffer = true;
  leftover_ind = 0;

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  for (; *cur_char != field_separator; ++cur_char) {
    // puts base either on buffer or leftover
    put_base_in_buffer();
    // Make sure we have memory space to process new base
    if (!write_to_buffer) {
      if (leftover_ind == stats.num_extra_bases) {
        stats.allocate_new_base(tile_ignore);
      }
    }

    // statistics updated base by base
    // use buffer
    if (write_to_buffer) {
      process_sequence_base_from_buffer(stats);
    }

    // use dynamic allocation
    else {
      process_sequence_base_from_leftover(stats);
    }
    // increase leftover position if not writing to buffer anymore
    if (!write_to_buffer) {
      leftover_ind++;
    }

    // either way increase read position
    ++read_pos;

    // if we reached the buffer size, stop
    if (read_pos == buffer_size) {
      write_to_buffer = false;
    }
  }

  // statistics summarized after the read
  postprocess_sequence_line(stats);
}

/*******************************************************/
/*************** QUALITY PROCESSING ********************/
/*******************************************************/
// Process quality value the fast way from buffer
inline void
StreamReader::process_quality_base_from_buffer(FastqStats &stats) {
  // Average quality in position
  stats.position_quality_count[
                               (read_pos << stats.kBitShiftQuality) | quality_value]++;

  // Tile processing
  if (!tile_ignore) {
    if ((stats.num_reads == next_tile_read) && tile_cur != 0) {
      stats.tile_position_quality[tile_cur][read_pos]
        += quality_value;
      stats.tile_position_count[tile_cur][read_pos]++;
    }
  }
}

// Slow version of function above
inline void
StreamReader::process_quality_base_from_leftover(FastqStats &stats) {
  // Average quality in position
  stats.long_position_quality_count[
                                    (leftover_ind << stats.kBitShiftQuality) | quality_value]++;

  // Tile processing
  if (!tile_ignore) {
    if ((stats.num_reads == next_tile_read) && tile_cur != 0) {
      stats.tile_position_quality[tile_cur][leftover_ind + stats.kNumBases]
        += quality_value;
      stats.tile_position_count[tile_cur][read_pos]++;
    }
  }
}

// Reads the quality line of each base.
inline void
StreamReader::read_quality_line(FastqStats &stats) {
  // reset quality counts
  read_pos = 0;
  cur_quality = 0;
  read_from_buffer = true;
  leftover_ind = 0;

  // For quality, we do not look for the separator, but rather for an explicit
  // newline or EOF in case the file does not end with newline or we are getting
  // decompressed strings from a stream
  for (; (*cur_char != line_separator) && !is_eof(); ++cur_char) {
    get_base_from_buffer();

    // Converts quality ascii to zero-based
    quality_value = *cur_char - stats.kBaseQuality;

    // Fast bases from buffer
    if (read_from_buffer) {
      process_quality_base_from_buffer(stats);
    }

    // Slow bases from dynamic allocation
    else {
      process_quality_base_from_leftover(stats);
    }

    // Sums quality value so we can bin the average at the end
    cur_quality += quality_value;

    if (!read_from_buffer) {
      ++leftover_ind;
    }

    // Flag to start reading and writing outside of buffer
    ++read_pos;
    if (read_pos == buffer_size) {
      read_from_buffer = false;
    }
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  if (do_quality_sequence) {
    stats.quality_count[cur_quality / read_pos]++;  // avg quality histogram
  }
}

/*******************************************************/
/*************** POST LINE PROCESSING ******************/
/*******************************************************/
/*************** THIS IS VERY SLOW ********************/
inline void
StreamReader::postprocess_fastq_record(FastqStats &stats) {
  if (do_overrepresented || do_duplication) {
    // if reads are >75pb, truncate to 50
    if (do_nogroup || read_pos <= stats.kDupReadMaxSize) {
      buffer[read_pos] = '\0';
    }
    else {
      buffer[stats.kDupReadTruncateSize] = '\0';
    }

    sequence_to_hash = string(buffer);

    // New sequence found
    if (stats.sequence_count.count(sequence_to_hash) == 0) {
      if (continue_storing_sequences) {
        stats.sequence_count.insert({{sequence_to_hash, 1}});
        stats.count_at_limit = stats.num_reads;
        ++stats.num_unique_seen;

        // if we reached the cutoff of 100k, stop storing
        if (stats.num_unique_seen == stats.kDupUniqueCutoff) {
          continue_storing_sequences = false;
        }
      }
    }
    else {
      stats.sequence_count[sequence_to_hash]++;
      if (continue_storing_sequences)
        stats.count_at_limit++;
    }
  }

  // counts tile if applicable
  if (!tile_ignore) {
    if (stats.num_reads == next_tile_read) {
      next_tile_read += num_reads_for_tile;
      if (tile_cur != 0) {
        stats.tile_count[tile_cur]++;
      }
    }
  }

  // I counted kmers here so register that I did so
  if (do_adapter) {
    if (stats.num_reads == next_kmer_read) {
      next_kmer_read += num_reads_for_kmer;
      stats.num_reads_kmer++;
    }
  }
}

/*******************************************************/
/*************** READ FASTQ RECORD *********************/
/*******************************************************/

// Set fastq separator as \n
FastqReader::FastqReader(FalcoConfig &_config,
                         const size_t _buffer_size) :
  StreamReader(_config, _buffer_size, '\n', '\n') {
}

bool inline
FastqReader::is_eof() {
  return cur_char == (last + 1);
}

// Load uncompressed fastq through memory map
void
FastqReader::load() {
  // uncompressed fastq = memorymap
  int fd = open(filename.c_str(), O_RDONLY, 0);
  if (fd == -1)
    throw runtime_error("failed to open fastq file: " + filename);

  // get the file size
  fstat(fd, &st);

  // execute mmap
  mmap_data = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (mmap_data == MAP_FAILED)
    throw runtime_error("failed to mmap fastq file: " + filename);

  // Initialize position pointer
  cur_char = static_cast<char*>(mmap_data);
  last = cur_char + st.st_size - 1;
}

// Parses the particular fastq format
bool inline
FastqReader::operator>>(FastqStats &stats) {
  read_tile_line(stats);
  skip_separator();

  read_sequence_line(stats);
  skip_separator();

  read_fast_forward_line();
  skip_separator();

  read_quality_line(stats);
  skip_separator();

  postprocess_fastq_record(stats);

  // Successful read, increment number in stats
  stats.num_reads++;

  // Returns if file should keep being checked
  return !is_eof();
}

FastqReader::~FastqReader()  {
  munmap(mmap_data, st.st_size);
}

#ifdef USE_ZLIB
/*******************************************************/
/*************** READ FASTQ GZ RCORD *******************/
/*******************************************************/
// the gz fastq constructor is the same as the fastq
GzFastqReader::GzFastqReader(FalcoConfig &_config,
                             const size_t _buffer_size) :
  FastqReader(_config, _buffer_size) {
}

// Load fastq with zlib
void
GzFastqReader::load() {
  fileobj = gzopen(filename.c_str(), "r");
  if (fileobj == Z_NULL)
    throw runtime_error("Cannot open gzip file : " + filename);
  cur_char = new char[1];
  ++cur_char;
}

// straightforward
bool inline
GzFastqReader::is_eof() {
  return gzeof(fileobj);
}

GzFastqReader::~GzFastqReader() {
  gzclose_r(fileobj);
}

// Parses fastq gz by reading line by line into the gzbuf
bool inline
GzFastqReader::operator >>(FastqStats &stats) {
  cur_char = gzgets(fileobj, gzbuf, kChunkSize);

  // need to check here if we did not hit eof
  if (is_eof()) {
    return false;
  }

  read_tile_line(stats);
  skip_separator();

  cur_char = gzgets(fileobj, gzbuf, kChunkSize);
  read_sequence_line(stats);
  skip_separator();

  cur_char = gzgets(fileobj, gzbuf, kChunkSize);
  read_fast_forward_line();
  skip_separator();

  cur_char = gzgets(fileobj, gzbuf, kChunkSize);
  read_quality_line(stats);
  skip_separator();

  postprocess_fastq_record(stats);

  // Successful read, increment number in stats
  stats.num_reads++;

  // Returns if file should keep being checked
  return !is_eof();
}
#endif

/*******************************************************/
/*************** READ SAM RECORD ***********************/
/*******************************************************/
// set sam separator as tab
SamReader::SamReader(FalcoConfig &_config,
                     const size_t _buffer_size) :
  StreamReader(_config, _buffer_size, '\t', '\n') {}

void
SamReader::load() {
  // uncompressed fastq = memorymap
  int fd = open(filename.c_str(), O_RDONLY, 0);
  if (fd == -1)
    throw runtime_error("failed to open fastq file: " + filename);

  // get the file size
  fstat(fd, &st);

  // execute mmap
  mmap_data = mmap(NULL, st.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (mmap_data == MAP_FAILED)
    throw runtime_error("failed to mmap fastq file: " + filename);

  // Initialize position pointer
  cur_char = static_cast<char*>(mmap_data);
  last = cur_char + st.st_size - 1;

  // Skip sam header
  while (*cur_char == '@') {
    for (; *cur_char != line_separator; ++cur_char) {}
    ++cur_char;
  }
}

bool inline
SamReader::is_eof() {
  return (cur_char == last + 1);
}

bool inline
SamReader::operator >> (FastqStats &stats) {
  read_tile_line(stats);
  skip_separator();
  for (size_t i = 0; i < 8; ++i) {
    read_fast_forward_line();
    skip_separator();
  }
  read_sequence_line(stats);
  read_quality_line(stats);

  // skips all tags after quality until newline
  while (*cur_char != line_separator) {
    read_fast_forward_line();
    skip_separator();
  }

  // skip \n
  ++cur_char;
  postprocess_fastq_record(stats);
  stats.num_reads++;

  // Returns if file should keep being checked
  return !is_eof();
}

SamReader::~SamReader() {
  munmap(mmap_data, st.st_size);
}

#ifdef USE_HTS
/*******************************************************/
/*************** READ BAM RECORD ***********************/
/*******************************************************/

// set sam separator as tab
BamReader::BamReader(FalcoConfig &_config, const size_t _buffer_size) :
  StreamReader(_config, _buffer_size, '\t', '\n') {
  rd_ret = 0;
}

void
BamReader::load() {
  if (!(hts = hts_open(filename.c_str(), "r")))
    throw runtime_error("cannot load bam file : " + filename);

  if (!(hdr = sam_hdr_read(hts)))
    throw runtime_error("failed to read header from file: " + filename);

  if (!(b = bam_init1()))
    throw runtime_error("failed to read record from file: " + filename);
}

// We will check eof on the >> operator
bool
BamReader::is_eof() {
  return (cur_char == last - 1);
}

bool inline
BamReader::operator>>(FastqStats &stats) {
  if ((rd_ret = sam_read1(hts, hdr, b)) >= 0) {
    fmt_ret = 0;
    if ((fmt_ret = sam_format1(hdr, b, &hts->line)) > 0) {
      // define char* values for processing lines char by char
      cur_char = hts->line.s;
      last = cur_char + strlen(hts->line.s) - 1;

      // Now read it as regular sam
      read_tile_line(stats);
      skip_separator();
      for (size_t i = 0; i < 8; ++i) {
        read_fast_forward_line();
        skip_separator();
      }

      read_sequence_line(stats);
      read_quality_line(stats);

      postprocess_fastq_record(stats);
      stats.num_reads++;
      return true;
    }
    else throw runtime_error("failed reading record from: " + filename);

    return false;
  }
  // If I could not read another line it means it's eof
  return false;
}

BamReader::~BamReader() {
  if (hdr) {
    bam_hdr_destroy(hdr);
    hdr = 0;
  }
  if (b) {
    bam_destroy1(b);
    b = 0;
  }
  if (hts) {
    hts_close(hts);
    hts = 0;
  }
}
#endif
