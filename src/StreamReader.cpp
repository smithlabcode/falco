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

#ifdef USE_HTS
#include "bamxx/bamxx.hpp"
#endif

#include <algorithm>
#include <cstring>
#include <fstream>
#include <string>
#include <vector>

// this is faster than std::min
template <class T>
inline T
min8(const T a, const T b) {
  return (a < b) ? a : b;
}
/****************************************************/
/***************** STREAMREADER *********************/
/****************************************************/
std::size_t
get_tile_split_position(const FalcoConfig &config) {
  const std::string &filename = config.filename;
  // Count colons to know the formatting pattern
  const infile_type_t infile_type = config.infile_type;
  std::string line;
  if (infile_type == infile_type_t::sam) {
    std::ifstream in(filename);
    if (!in)
      throw std::runtime_error("cannot load sam file : " + filename);
    while (std::getline(in, line) && std::size(line) > 0 && line[0] == '@')
      ;
    const auto first_tab = std::find(std::cbegin(line), std::cend(line), '\t');
    line.resize(std::distance(std::cbegin(line), first_tab));
  }
#ifdef USE_HTS
  if (infile_type == infile_type_t::bam) {
    bamxx::bam_in hts(filename);
    if (!hts)
      throw std::runtime_error("failed to open BAM file: " + filename);
    bamxx::bam_header hdr(hts);
    if (!hdr)
      throw std::runtime_error("failed to read BAM header: " + filename);
    bamxx::bam_rec aln;
    if (!hts.read(hdr, aln))
      throw std::runtime_error("failed to read BAM record from: " + filename);
    line = std::string{bam_get_qname(aln.b)};
  }
#endif
  if (infile_type == infile_type_t::fastq_gz) {
    gzFile in = gzopen(std::data(filename), "rb");
    if (!in)
      throw std::runtime_error("problem reading input file: " + filename);
    // Read the first line of the file
    static constexpr auto buffer_size = 1024u;
    line.resize(buffer_size, '\0');
    gzgets(in, std::data(line), buffer_size);
    gzclose(in);
  }
  if (infile_type == infile_type_t::fastq) {
    std::ifstream in(filename);
    if (!in)
      throw std::runtime_error("problem reading input file: " + filename);
    std::getline(in, line);
  }

  const auto num_colon = std::count(std::cbegin(line), std::cend(line), ':');

  // Copied from fastqc
  if (num_colon >= 6)
    return 4;
  if (num_colon >= 4)
    return 2;
  return 0;  // no tile information on read name
}

// function to turn a vector into array for adapter hashes and faster lookup
std::array<std::size_t, Constants::max_adapters>
make_adapters(const std::vector<std::size_t> &adapter_hashes) {
  static constexpr auto error_message =
    "Number of adapters is larger than 128, which hinders visualziation and "
    "speed of falco. Please keep it to under 128";
  if (std::size(adapter_hashes) > Constants::max_adapters)
    throw std::runtime_error(error_message);
  std::array<std::size_t, Constants::max_adapters> ans;
  std::copy(std::cbegin(adapter_hashes), std::cend(adapter_hashes),
            std::begin(ans));
  return ans;
}

StreamReader::StreamReader(FalcoConfig &config, const std::size_t buffer_size,
                           const char field_separator,
                           const char line_separator) :
  // must pass the config skips as const to read them fast
  do_sequence_hash(config.do_duplication || config.do_overrepresented),
  do_kmer(config.do_kmer), do_adapter(config.do_adapter),
  do_adapter_optimized(config.do_adapter_optimized),
  do_sliding_window(do_adapter_optimized || do_kmer),
  do_n_content(config.do_n_content), do_quality_base(config.do_quality_base),
  do_sequence(config.do_sequence), do_gc_sequence(config.do_gc_sequence),
  do_quality_sequence(config.do_quality_sequence), do_tile(config.do_tile),
  do_sequence_length(config.do_sequence_length),

  // Here are the usual stream reader configs
  field_separator(field_separator), line_separator(line_separator),
  buffer_size(buffer_size), read_step(config.read_step),
  tile_split_point(get_tile_split_position(config)),
  tile_ignore(!do_tile || tile_split_point == 0),

  // Here are the const adapters
  do_adapters_slow(config.do_adapter && !config.do_adapter_optimized),
  adapter_seqs(config.adapter_seqs),

  num_adapters(std::size(config.adapter_hashes)),
  adapter_size(config.adapter_size),
  // for case size == 32 expr (1ull << 64) -1 gives 0.
  // We need to set mask as all 64 bits 1 => use SIZE_MAX in this case
  adapter_mask(adapter_size == 32 ? SIZE_MAX
                                  : (1ull << (2 * adapter_size)) - 1),
  adapters(make_adapters(config.adapter_hashes)), filename(config.filename) {

  // Allocates buffer to temporarily store reads
  buffer = new char[buffer_size + 1];  // +1 for the \0
  buffer[buffer_size] = '\0';

  // duplication init
  continue_storing_sequences = true;

  // Tile init
  tile_cur = 0;

  // keep track of which reads to do tile
  next_read = 0;
  next_tile_read = 0;
  next_kmer_read = 0;
  do_tile_read = true;

  // GS: test
  leftover_ind = 0;
}

// Makes sure that any subclass deletes the buffer
StreamReader::~StreamReader() { delete[] buffer; }

/*******************************************************/
/*************** BUFFER MANAGEMENT *********************/
/*******************************************************/
// puts base either on buffer or leftover
void
StreamReader::put_base_in_buffer() {
  base_from_buffer = *cur_char;
  if (still_in_buffer) {
    buffer[read_pos] = base_from_buffer;
    return;
  }
  if (leftover_ind > std::size(leftover_buffer))
    leftover_buffer.resize(leftover_ind);
  leftover_buffer[leftover_ind] = base_from_buffer;
}

// Gets base from either buffer or leftover
void
StreamReader::get_base_from_buffer() {
  base_from_buffer =
    still_in_buffer ? buffer[read_pos] : leftover_buffer[leftover_ind];
}

/*******************************************************/
/*************** FAST FOWARD ***************************/
/*******************************************************/

// Keeps going forward while the current character is a separator
void
StreamReader::skip_separator() {
  while (*cur_char == field_separator)
    ++cur_char;
}

// Skips lines that are not relevant
void
StreamReader::read_fast_forward_line() {
  while (*cur_char != field_separator)
    ++cur_char;
}

void
StreamReader::read_fast_forward_line_eof() {
  while (*cur_char != field_separator && *cur_char != line_separator &&
         !is_eof())
    ++cur_char;
}

/*******************************************************/
/*************** TILE PROCESSING ***********************/
/*******************************************************/
// Parse the comment

void
StreamReader::get_tile_value() {
  tile_cur = 0;
  std::size_t num_colon = 0;
  for (; *cur_char != field_separator; ++cur_char) {
    num_colon += (*cur_char == ':');
    if (num_colon == tile_split_point) {
      ++cur_char;  // skip colon
      // parse till next colon or \n
      for (; *cur_char != ':' && *cur_char != field_separator; ++cur_char)
        tile_cur = tile_cur * 10 + (*cur_char - '0');
      // now fast forward until last \n
      while (*cur_char != field_separator)
        ++cur_char;
      return;
    }
  }
}

// Gets the tile from the sequence name (if applicable)
void
StreamReader::read_tile_line(FastqStats &stats) {
  do_tile_read = (do_read && stats.num_reads == next_tile_read);
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

  get_tile_value();
  // allocate vector for tile if it doesn't exist
  if (stats.tile_position_quality.find(tile_cur) ==
      std::cend(stats.tile_position_quality)) {
    stats.tile_position_quality[tile_cur] =
      std::vector<double>(stats.max_read_length, 0.0);
    // stats.tile_position_quality.find(tile_cur)->second[0] = 0;
    stats.tile_position_count[tile_cur] =
      std::vector<std::size_t>(stats.max_read_length, 0);
  }
}

/*******************************************************/
/*************** SEQUENCE PROCESSING *******************/
/*******************************************************/

// This is probably the most important function for speed, so it must be really
// optimized at all times
void
StreamReader::process_sequence_base_from_buffer(FastqStats &stats) {
  // count Ns even if asked not to report them
  if (base_from_buffer == 'N') {
    ++stats.n_base_count[read_pos];
    num_bases_after_n = 1;  // start over the current kmer
    return;
  }

  const auto two_bit = actg_to_2bit(base_from_buffer);

  // ATGC bases, increments basic statistic counts
  cur_gc_count += (two_bit & 1);
  ++stats.base_count[(read_pos << Constants::bit_shift_base) | two_bit];

  if (do_sliding_window) {  // Update k-mer sequence
    cur_kmer = (cur_kmer << Constants::bit_shift_base) | two_bit;
    // registers k-mer if seen at least k nucleotides since the last n
    if (do_kmer && do_kmer_read && num_bases_after_n >= Constants::kmer_size) {
      ++stats.kmer_count[(read_pos << Constants::bit_shift_kmer) |
                         (cur_kmer & Constants::kmer_mask)];
      ++stats.pos_kmer_count[read_pos];
    }

    // GS: slow, need to use fsm
    if (do_adapter_optimized && num_bases_after_n == adapter_size) {
      cur_kmer &= adapter_mask;
      for (std::size_t i = 0; i < num_adapters; ++i)
        if (cur_kmer == adapters[i] && !adapters_found[i]) {
          const auto shifted = read_pos << Constants::bit_shift_adapter;
          ++stats.pos_adapter_count[shifted | i];
          adapters_found[i] = true;
        }
    }
    num_bases_after_n += (num_bases_after_n != adapter_size);
  }
}

// slower version of process_sequence_base_from_buffer that dynamically
// allocates if base position is not already cached
void
StreamReader::process_sequence_base_from_leftover(FastqStats &stats) {
  if (base_from_buffer == 'N') {
    ++stats.long_n_base_count[leftover_ind];
    num_bases_after_n = 1;  // start over the current kmer
    return;
  }
  // ATGC bases
  // increments basic statistic counts
  const auto two_bit = actg_to_2bit(base_from_buffer);
  cur_gc_count += (two_bit & 1);
  const auto idx = (leftover_ind << Constants::bit_shift_base) | two_bit;
  ++stats.long_base_count[idx];
  // WE WILL NOT DO KMER STATS OUTSIDE OF BUFFER
}

// Gets statistics after reading the entire sequence line
void
StreamReader::postprocess_sequence_line(FastqStats &stats) {
  // Updates basic statistics total GC
  stats.total_gc += cur_gc_count;

  // read length frequency histogram
  if (do_sequence_length) {
    if (still_in_buffer) {
      stats.empty_reads += (read_pos == 0);
      if (read_pos != 0)
        ++stats.read_length_freq[read_pos - 1];
    }
    else
      ++stats.long_read_length_freq[leftover_ind - 1];
  }

  // Updates maximum read length if applicable
  stats.max_read_length =
    read_pos > stats.max_read_length ? read_pos : stats.max_read_length;
  stats.total_bases += read_pos;

  // FastQC's gc model summarized, if requested
  if (do_gc_sequence && read_pos != 0) {
    // If we haven't passed the short base threshold, we use the cached models
    if (still_in_buffer) {
      // if we haven't passed the truncation point, use the current values,
      // otherwise we have truncated previously
      if (next_truncation == 100) {
        truncated_length = read_pos;
        truncated_gc_count = cur_gc_count;
      }
      for (const auto &v :
           stats.gc_models[truncated_length].models[truncated_gc_count]) {
        stats.gc_count[v.percent] += v.increment;
      }

      // if the read length is too large, we just use the discrete percentage
    }
    else {
      ++stats.gc_count[100 * cur_gc_count / read_pos];
    }
  }
}

// Reads the line that has the biological sequence
void
StreamReader::read_sequence_line(FastqStats &stats) {
  if (!do_read) {
    read_fast_forward_line();
    return;
  }
  // restart line counters
  read_pos = 0;
  cur_gc_count = 0;
  truncated_gc_count = 0;
  num_bases_after_n = 1;
  still_in_buffer = true;
  next_truncation = 100;
  do_kmer_read = (stats.num_reads == next_kmer_read);
  adapters_found.reset();

  if (do_adapters_slow) {
    const std::string seq_line_str = cur_char;
    for (std::size_t i = 0; i != num_adapters; ++i) {
      const std::size_t adapt_index = seq_line_str.find(adapter_seqs[i]);
      if (adapt_index < stats.SHORT_READ_THRESHOLD) {
        const auto x = ((adapt_index + std::size(adapter_seqs[i]) - 1)
                        << Constants::bit_shift_adapter);
        ++stats.pos_adapter_count[x | i];
      }
    }
  }

  /*********************************************************/
  /********** THIS LOOP MUST BE ALWAYS OPTIMIZED ***********/
  /*********************************************************/
  for (; *cur_char != field_separator; ++cur_char, ++read_pos) {
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

    // Truncate GC counts to multiples of 100
    if (do_gc_sequence && read_pos == next_truncation) {
      truncated_gc_count = cur_gc_count;
      truncated_length = read_pos;
      next_truncation += 100;
    }
  }

  // statistics summarized after the read
  postprocess_sequence_line(stats);
}

/*******************************************************/
/*************** QUALITY PROCESSING ********************/
/*******************************************************/
// Process quality value the fast way from buffer
void
StreamReader::process_quality_base_from_buffer(FastqStats &stats) {
  // Average quality in position
  ++stats.position_quality_count[(read_pos << Constants::bit_shift_quality) |
                                 quality_value];

  // Tile processing
  if (!tile_ignore && do_tile_read && tile_cur != 0) {
    // allocate more base space if necessary
    if (std::size(stats.tile_position_quality[tile_cur]) == read_pos) {
      stats.tile_position_quality[tile_cur].push_back(0.0);
      stats.tile_position_count[tile_cur].push_back(0);
    }
    stats.tile_position_quality[tile_cur][read_pos] += quality_value;
    ++stats.tile_position_count[tile_cur][read_pos];
  }
}

// Slow version of function above
void
StreamReader::process_quality_base_from_leftover(FastqStats &stats) {
  // Average quality in position
  ++stats.long_position_quality_count[(leftover_ind
                                       << Constants::bit_shift_quality) |
                                      quality_value];

  // Tile processing
  if (!tile_ignore) {
    if (do_tile_read && tile_cur != 0) {
      stats.tile_position_quality[tile_cur][read_pos] += quality_value;
      ++stats.tile_position_count[tile_cur][read_pos];
    }
  }
}

// Reads the quality line of each base.
void
StreamReader::read_quality_line(FastqStats &stats) {
  // MN: Below, the second condition tests whether the quality score in sam
  // file only consists of '*', which indicates that no quality score is
  // available
  if (!do_read || (*cur_char == '*' && (*(cur_char + 1) == field_separator ||
                                        *(cur_char + 1) == field_separator))) {
    read_fast_forward_line_eof();
    return;
  }
  // reset quality counts
  read_pos = 0;
  cur_quality = 0;
  still_in_buffer = true;

  // For quality, we do not look for the separator, but rather for an explicit
  // newline or EOF in case the file does not end with newline or we are getting
  // decompressed strings from a stream
  while (*cur_char != field_separator && *cur_char != line_separator &&
         !is_eof()) {
    if (read_pos == buffer_size) {
      still_in_buffer = false;
      leftover_ind = 0;
    }

    get_base_from_buffer();

    // update lowest quality
    stats.lowest_char =
      min8(stats.lowest_char,
           *cur_char == 9 ? std::numeric_limits<char>::max() : *cur_char);

    // Converts quality ascii to zero-based
    quality_value = *cur_char - Constants::quality_zero;

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

    ++cur_char;  // increment loop
  }

  // Average quality approximated to the nearest integer. Used to make a
  // histogram in the end of the summary.
  if (read_pos != 0)
    ++stats.quality_count[cur_quality / read_pos];  // avg quality histogram
}

/*******************************************************/
/*************** POST LINE PROCESSING ******************/
/*******************************************************/
/*************** THIS IS VERY SLOW ********************/
// if reads are >75pb, truncate to 50 akin to FastQC
inline std::size_t
get_truncate_point(const std::size_t read_pos) {
  return (read_pos <= Constants::unique_reads_max_length)
           ? read_pos
           : Constants::unique_reads_truncate;
}

void
StreamReader::postprocess_fastq_record(FastqStats &stats) {
  if (do_sequence_hash) {
    buffer[get_truncate_point(read_pos)] = '\0';
    sequence_to_hash = std::string(buffer);
    // New sequence found
    if (stats.sequence_count.count(sequence_to_hash) == 0) {
      if (continue_storing_sequences) {
        stats.sequence_count.insert({{sequence_to_hash, 1}});
        stats.count_at_limit = stats.num_reads;
        ++stats.num_unique_seen;

        // if we reached the cutoff of 100k, stop storing
        if (stats.num_unique_seen == Constants::unique_reads_stop_counting)
          continue_storing_sequences = false;
      }
    }
    else {
      ++stats.sequence_count[sequence_to_hash];
      stats.count_at_limit += continue_storing_sequences;
    }
  }
  // counts tile if applicable
  if (!tile_ignore) {
    if (do_tile_read) {
      next_tile_read += num_reads_for_tile;
    }
  }
  next_kmer_read += do_kmer_read * num_reads_for_kmer;
}

inline bool
StreamReader::check_bytes_read(const std::size_t read_num) {
  return (read_num & check_bytes_read_mask) == 0;
}

// read FASTQ record
char
get_line_separator(const std::string &filename) {
  auto fp = std::fopen(std::data(filename), "r");
  if (fp == nullptr)
    throw std::runtime_error("bad input file: " + filename);
  while (!std::feof(fp)) {
    const char c = std::fgetc(fp);
    if (c == '\n' || c == '\r') {
      std::fclose(fp);
      return c;
    }
  }
  std::fclose(fp);
  return '\n';
}

// Set fastq field_separator as line_separator
FastqReader::FastqReader(FalcoConfig &config, const std::size_t buffer_size) :
  StreamReader(config, buffer_size, get_line_separator(config.filename),
               get_line_separator(config.filename)) {
  file_buffer.resize(RESERVE_SIZE);
}

std::size_t
get_file_size(const std::string &filename) {
  auto fp = std::fopen(std::data(filename), "r");
  if (fp == nullptr)
    throw std::runtime_error("bad input file: " + filename);
  std::fseek(fp, 0L, SEEK_END);
  const auto ret = std::ftell(fp);
  std::fclose(fp);
  return ret;
}

// Load fastq with zlib
[[nodiscard]] std::size_t
FastqReader::load() {
  fileobj = std::fopen(std::data(filename), "r");
  if (fileobj == NULL)
    throw std::runtime_error("Cannot open FASTQ file : " + filename);
  return get_file_size(filename);
}

// straightforward
inline bool
FastqReader::is_eof() {
  return std::feof(fileobj);
}

FastqReader::~FastqReader() { std::fclose(fileobj); }

// Parses fastq gz by reading line by line into the gzbuf
bool
FastqReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  cur_char = std::fgets(std::data(file_buffer), RESERVE_SIZE, fileobj);

  // need to check here if we did not hit eof
  if (is_eof())
    return false;

  do_read = (stats.num_reads == next_read);

  read_tile_line(stats);
  skip_separator();

  cur_char = std::fgets(std::data(file_buffer), RESERVE_SIZE, fileobj);
  read_sequence_line(stats);

  skip_separator();

  cur_char = std::fgets(std::data(file_buffer), RESERVE_SIZE, fileobj);
  read_fast_forward_line();
  skip_separator();

  cur_char = std::fgets(std::data(file_buffer), RESERVE_SIZE, fileobj);
  read_quality_line(stats);
  skip_separator();

  if (do_read)
    postprocess_fastq_record(stats);

  next_read += do_read * read_step;

  // Successful read, increment number in stats
  ++stats.num_reads;

  // Returns if file should keep being checked
  if (check_bytes_read(stats.num_reads))
    num_bytes_read = std::ftell(fileobj);

  return (!is_eof() && cur_char != 0);
}

/*******************************************************/
/*************** READ FASTQ GZ RCORD *******************/
/*******************************************************/
// the gz fastq constructor is the same as the fastq
GzFastqReader::GzFastqReader(FalcoConfig &config,
                             const std::size_t buffer_size) :
  StreamReader(config, buffer_size, '\n', '\n') {
  gzbuf.resize(RESERVE_SIZE);
}

// Load fastq with zlib
[[nodiscard]] std::size_t
GzFastqReader::load() {
  fileobj = gzopen(std::data(filename), "r");
  if (fileobj == Z_NULL)
    throw std::runtime_error("Cannot open gzip FASTQ file : " + filename);

  return get_file_size(filename);
}

// straightforward
inline bool
GzFastqReader::is_eof() {
  return gzeof(fileobj);
}

GzFastqReader::~GzFastqReader() { gzclose_r(fileobj); }

// Parses fastq gz by reading line by line into the gzbuf
bool
GzFastqReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  cur_char = gzgets(fileobj, std::data(gzbuf), RESERVE_SIZE);

  // need to check here if we did not hit eof
  if (is_eof()) {
    return false;
  }
  do_read = (stats.num_reads == next_read);

  read_tile_line(stats);
  skip_separator();

  cur_char = gzgets(fileobj, std::data(gzbuf), RESERVE_SIZE);
  read_sequence_line(stats);
  skip_separator();

  cur_char = gzgets(fileobj, std::data(gzbuf), RESERVE_SIZE);
  read_fast_forward_line();
  skip_separator();

  cur_char = gzgets(fileobj, std::data(gzbuf), RESERVE_SIZE);
  read_quality_line(stats);
  skip_separator();

  if (do_read)
    postprocess_fastq_record(stats);

  next_read += do_read * read_step;

  // Successful read, increment number in stats
  ++stats.num_reads;

  // Returns if file should keep being checked
  if (check_bytes_read(stats.num_reads))
    num_bytes_read = gzoffset(fileobj);
  return !is_eof();
}

/*******************************************************/
/*************** READ SAM RECORD ***********************/
/*******************************************************/
// set sam separator as tab
SamReader::SamReader(FalcoConfig &config, const std::size_t buffer_size) :
  StreamReader(config, buffer_size, '\t', get_line_separator(config.filename)) {
  file_buffer.resize(RESERVE_SIZE);
}

[[nodiscard]] std::size_t
SamReader::load() {
  fileobj = std::fopen(std::data(filename), "r");
  if (fileobj == NULL)
    throw std::runtime_error("Cannot open SAM file : " + filename);

  // skip sam header
  char first_char_in_line;
  while (!is_eof() && ((first_char_in_line = std::fgetc(fileobj)) == '@')) {
    std::ungetc(first_char_in_line, fileobj);
    cur_char = std::fgets(std::data(file_buffer), RESERVE_SIZE, fileobj);
  }
  return get_file_size(filename);
}

inline bool
SamReader::is_eof() {
  return std::feof(fileobj);
}

bool
SamReader::read_entry(FastqStats &stats, std::size_t &num_bytes_read) {
  cur_char = std::fgets(std::data(file_buffer), RESERVE_SIZE, fileobj);

  if (is_eof())
    return false;
  do_read = (stats.num_reads == next_read);

  read_tile_line(stats);
  skip_separator();

  for (std::size_t i = 0; i < 8; ++i) {
    read_fast_forward_line();
    skip_separator();
  }

  // field 10
  read_sequence_line(stats);
  skip_separator();

  // field 11
  read_quality_line(stats);
  if (do_read)
    postprocess_fastq_record(stats);

  next_read += do_read * read_step;
  ++stats.num_reads;

  // skips all tags after quality until newline
  while (*cur_char != line_separator && !is_eof()) {
    if (check_bytes_read(stats.num_reads))
      num_bytes_read = std::ftell(fileobj);
    ++cur_char;
  }

  // Returns if file should keep being checked
  return !is_eof() && cur_char != 0;
}

SamReader::~SamReader() { std::fclose(fileobj); }
