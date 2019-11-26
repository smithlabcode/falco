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
using std::array;
/****************************************************/
/***************** STREAMREADER *********************/
/****************************************************/

// function to turn a vector into array for adapter hashes and fast lookup
array<size_t, FastqStats::max_adapters> 
make_adapters(const vector<size_t> &adapter_hashes) {
  if (adapter_hashes.size() > FastqStats::max_adapters)
    throw runtime_error("Number of adapters is larger than 128, which hinders "
                        "visualziation and speed of falco. Please keep it to "
                        "under 128");

  array<size_t, FastqStats::max_adapters> ans;
  for (size_t i = 0; i < adapter_hashes.size(); ++i)
    ans[i] = adapter_hashes[i];

  return ans;
}

StreamReader::StreamReader(FalcoConfig &config,
                           const size_t _buffer_size,
                           const char _field_separator,
                           const char _line_separator) :
  // I have to pass the config skips as const to read them fast
  do_sequence_hash(config.do_duplication || config.do_overrepresented),
  do_kmer(config.do_kmer),
  do_adapter(config.do_adapter),
  do_sliding_window(do_adapter || do_kmer),
  do_n_content(config.do_n_content),
  do_quality_base(config.do_quality_base),
  do_sequence(config.do_sequence),
  do_gc_sequence(config.do_gc_sequence),
  do_quality_sequence(config.do_quality_sequence),
  do_tile(config.do_tile),
  do_sequence_length(config.do_sequence_length),

  // Here are the usual stream reader configs
  buffer_size(_buffer_size),
  field_separator(_field_separator),
  line_separator(_line_separator),

  // Here are the const adapters
  num_adapters(config.adapter_hashes.size()),
  adapters(make_adapters(config.adapter_hashes))
  {

  // Allocates buffer to temporarily store reads
  buffer = new char[buffer_size + 1];
  buffer[buffer_size] = '\0';

  // duplication init
  continue_storing_sequences = true;

  // Tile init
  tile_ignore = !do_tile;  // early ignore tile if asked to skip it
  tile_cur = 0;
  tile_split_point = 0;

  // keep track of which reads to do tile
  next_tile_read = 0;
  do_tile_read = true;

  // Subclasses will use this to deflate if necessary
  filename = config.filename;

  // GS: test
  leftover_ind = 0;
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
  if (still_in_buffer) {
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
  if (still_in_buffer) {
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
    num_colon += (*cur_char == ':');
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
    num_colon += (*cur_char == ':');
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

  do_tile_read = (stats.num_reads == next_tile_read);
  if (!do_tile_read) {
    read_fast_forward_line();
    return;
  }
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
    stats.base_count[
      (read_pos << stats.kBitShiftNucleotide) | base_ind]++;

    if (do_sliding_window) {
      // Update k-mer sequence
      cur_kmer = ((cur_kmer << stats.kBitShiftNucleotide) | base_ind);

      // registers k-mer if seen at least k nucleotides since the last n
      if (do_kmer && (num_bases_after_n == stats.kmer_size)) {

          stats.kmer_count[(read_pos << stats.kBitShiftKmer)
                           | (cur_kmer & stats.kmer_mask)]++;
          stats.pos_kmer_count[read_pos]++;
      }

      // GS: slow, need to use fsm
      if (do_adapter && (num_bases_after_n == stats.adapter_size)) {
        cur_kmer &= stats.adapter_mask;
        for (i = 0; i < num_adapters; ++i) {
          if (cur_kmer == adapters[i]) {
            stats.pos_adapter_count[(read_pos << stats.kBitShiftAdapter) | i]++;
          }
        }
      }

      num_bases_after_n += (num_bases_after_n != stats.adapter_size);
    }
  }
}

// slower version of process_sequence_base_from_buffer that dynamically
// allocates if base position is not already cached
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
  // Updates basic statistics total GC
  stats.total_gc += cur_gc_count;

  // read length frequency histogram
  if (do_sequence_length) {
    if (still_in_buffer) {
      stats.read_length_freq[read_pos - 1]++;
    } else {
      stats.long_read_length_freq[leftover_ind - 1]++;
    }
  }

  // Updates maximum read length if applicable
  if (read_pos > stats.max_read_length) {
    stats.max_read_length = read_pos;
  }

  // FastQC's gc model summarized, if requested
  if (do_gc_sequence) {
    // If we haven't passed the short base threshold, we use the cached models
    if (still_in_buffer) {
      // if we haven't passed the truncation point, use the current values,
      // otherwise we have truncated previously
      if (next_truncation == 100) {
        truncated_length = read_pos;
        truncated_gc_count = cur_gc_count;
      }
      for (auto v :
          stats.gc_models[truncated_length].models[truncated_gc_count]) {
        stats.gc_count[v.percent] += v.increment;
      }

    // if the read length is too large, we just use the discrete percentage
    } else {
      stats.gc_count[100 * cur_gc_count / read_pos]++;
    }
  }
}

// Reads the line that has the biological sequence
inline void
StreamReader::read_sequence_line(FastqStats &stats) {
  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  truncated_gc_count = 0;
  num_bases_after_n = 1;
  still_in_buffer = true;
  next_truncation = 100;

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  for (; *cur_char != field_separator; ++cur_char) {
    // if we reached the buffer size, stop using it and start using leftover
    if (read_pos == buffer_size) {
      still_in_buffer = false;
      leftover_ind = 0;
    }

    // Make sure we have memory space to process new base
    if (!still_in_buffer) {
      if (leftover_ind == stats.num_extra_bases) {
        stats.allocate_new_base(tile_ignore);
      }
    }

    // puts base either on buffer or leftover
    put_base_in_buffer();

    // statistics updated base by base
    // use buffer
    if (still_in_buffer) {
      process_sequence_base_from_buffer(stats);
    }

    // use dynamic allocation
    else {
      process_sequence_base_from_leftover(stats);

      // Increase leftover pos if no longer in buffer
      ++leftover_ind;
    }

    // either way increase read position
    ++read_pos;


    // Truncate GC counts to multiples of 100
    if (do_gc_sequence) {
      if (read_pos == next_truncation) {
        truncated_gc_count = cur_gc_count;
        truncated_length= read_pos;
        next_truncation += 100;
      }
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
    (read_pos << stats.kBitShiftQuality) | quality_value
  ]++;

  // Tile processing
  if (!tile_ignore) {
    if (do_tile_read && tile_cur != 0) {
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
    if (do_tile_read && tile_cur != 0) {
      stats.tile_position_quality[tile_cur][read_pos]
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
  still_in_buffer = true;

  // For quality, we do not look for the separator, but rather for an explicit
  // newline or EOF in case the file does not end with newline or we are getting
  // decompressed strings from a stream
  for (; (*cur_char != line_separator) && !is_eof(); ++cur_char) {
    if (read_pos == buffer_size) {
      still_in_buffer = false;
      leftover_ind = 0;
    }

    get_base_from_buffer();

    // Converts quality ascii to zero-based
    quality_value = *cur_char - stats.kBaseQuality;

    // Fast bases from buffer
    if (still_in_buffer) {
      process_quality_base_from_buffer(stats);
    }

    // Slow bases from dynamic allocation
    else {
      process_quality_base_from_leftover(stats);
      ++leftover_ind;
    }

    // Sums quality value so we can bin the average at the end
    cur_quality += quality_value;

    // Flag to start reading and writing outside of buffer
    ++read_pos;
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  stats.quality_count[cur_quality / read_pos]++;  // avg quality histogram
}

/*******************************************************/
/*************** POST LINE PROCESSING ******************/
/*******************************************************/
/*************** THIS IS VERY SLOW ********************/
inline void
StreamReader::postprocess_fastq_record(FastqStats &stats) {
  if (do_sequence_hash) {
    // if reads are >75pb, truncate to 50
    if (read_pos <= stats.kDupReadMaxSize) {
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
      stats.count_at_limit += continue_storing_sequences;
    }
  }
  // counts tile if applicable
  if (!tile_ignore) {
    if (do_tile_read) {
      next_tile_read += num_reads_for_tile;
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
