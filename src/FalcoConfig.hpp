/* Copyright (C) 2019-2026 Guilherme De Sena Brandine and
 *                         Andrew D. Smith
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

#ifndef FALCO_CONFIG_HPP_
#define FALCO_CONFIG_HPP_

#include "aux.hpp"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

enum class infile_type_t : std::uint8_t {
  sam = 0,
  bam = 1,
  fastq = 2,
  fastq_gz = 3,
};

[[nodiscard]] inline std::string
to_string(const infile_type_t x) {
  switch (x) {
  case infile_type_t::sam:
    return "SAM";
  case infile_type_t::bam:
    return "BAM";
  case infile_type_t::fastq:
    return "FASTQ";
  case infile_type_t::fastq_gz:
    return "compressed FASTQ";
  default:
    return "unknown type";
  }
}

// config from options, constants, magic numbers, etc
struct FalcoConfig {
  FalcoConfig(const int argc, char *argv[]);

  // FASTQC OPTION PARSER
  bool casava{};            // files from raw casava output
  bool nanopore{};          // fast5 format
  bool nofilter{};          // if running with --casava flag
  bool extract{};           // if set the zipped file will be uncompressed
  bool nogroup{};           // disable grouping of bases for reads >50bp
  bool compressed{};        // whether or not to inflate file
  bool quiet{};             // suppress all progress output to terminal
  std::size_t read_step{};  // only process reads that are multiple of read_step
  std::uint32_t threads{1};       // number of threads to read files in parallel
  std::string call;               // the function call
  std::string format;             // force file format
  std::string contaminants_file;  // custom contaminants file
  std::string adapters_file;      // adapters file
  std::string limits_file;  // file with limits and options and custom analyses
  static const std::string html_template;  // the html for the template
  std::string tmpdir;  // dir for temp files when generating report images

  // Falco only
  bool progress{};  // report the progress bar

  // config on how to handle reads
  bool do_duplication{};
  bool do_kmer{};
  bool do_n_content{};
  bool do_overrepresented{};
  bool do_quality_base{};
  bool do_sequence{};
  bool do_gc_sequence{};
  bool do_quality_sequence{};
  bool do_tile{};
  bool do_adapter{};
  bool do_adapter_optimized{};
  bool do_sequence_length{};

  // Fastqc limits
  // These will become const bools in the stream reader
  std::unordered_map<std::string, std::unordered_map<std::string, double>>
    limits;
  // static const std::vector<std::string> values_to_check;
  static const std::unordered_set<std::string> values_to_check;

  // Contaminants
  // below: first = name, scond = seq
  std::vector<std::pair<std::string, std::string>> contaminants;

  /*************** ADAPTERS *********************/
  // Name (eg: Illumina Small RNA adapter)
  std::vector<std::string> adapter_names;

  // Actual string sequence (eg: ATTGCCACA)
  std::vector<std::string> adapter_seqs;

  // two-bit hash of the sequence above
  std::vector<std::size_t> adapter_hashes;

  std::size_t adapter_size{};
  std::size_t shortest_adapter_size{};

  // Additional information about the sample
  bool is_bisulfite{};
  bool is_reverse_complement{};

  infile_type_t infile_type{};

  std::string filename;
  std::string filename_stripped;

  void
  define_file_format();

  void
  read_limits();

  void
  read_adapters();

  void
  read_contaminants_file();

  void
  setup();
};

/*************************************************************
 ******************** ALL MAGIC NUMBERS **********************
 *************************************************************/
namespace Constants {
// log of a power of two, to use in bit shifting for fast index acces returns
// the log2 of a number if it is a power of two, or zero otherwise

// clang-format off
template<typename T>
constexpr T
log2exact(const T v) {
  return 63 -
    (v & 0x00000000FFFFFFFF ? 32 : 0) -
    (v & 0x0000FFFF0000FFFF ? 16 : 0) -
    (v & 0x00FF00FF00FF00FF ?  8 : 0) -
    (v & 0x0F0F0F0F0F0F0F0F ?  4 : 0) -
    (v & 0x3333333333333333 ?  2 : 0) -
    (v & 0x5555555555555555 ?  1 : 0);
}
// clang-format on

static constexpr std::size_t kmer_size = 7;
static constexpr std::size_t max_adapters = 128;

// number of bases for static allocation.
static constexpr std::size_t num_static_bases = 500;

// Value to subtract quality characters to get the actual quality value
static constexpr std::size_t quality_zero = 33;  // ascii for the lowest quality

// Smallest power of two that comprises all possible Illumina quality values.
// Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of is
// to avoid double pointer jumps and to get indices with bit shifts.
static constexpr std::size_t num_quality_values = 128;

// How many possible nucleotides (must be power of 2!)
static constexpr std::size_t num_nucleotides = 4;  // A=00,C=01,T=10,G=11

/************* DUPLICATION ESTIMATES *************/
// Number of unique sequences to see before stopping counting sequences
static constexpr std::size_t unique_reads_stop_counting = 1e5;

// Maximum read length to store the entire read in memory
static constexpr std::size_t unique_reads_max_length = 75;

// Prefix size to cut if read length exceeds the value above
static constexpr std::size_t unique_reads_truncate = 50;

// Bit shifts as instructions for the arrays
// for matrices that count stats per nucleotide
static constexpr std::size_t bit_shift_base = log2exact(num_nucleotides);

// for matrices that count stats for quality value
static constexpr std::size_t bit_shift_quality = log2exact(num_quality_values);

// bit shift for adapters, log(128)=7
static constexpr std::size_t bit_shift_adapter = log2exact(max_adapters);

// we shift 14 bits when reading a kmer, two bits per base
static constexpr std::size_t bit_shift_kmer = bit_shift_base * kmer_size;

// mask to get only the first 2*k bits of the sliding window
static constexpr std::size_t kmer_mask = (1ull << (bit_shift_kmer)) - 1;
};  // namespace Constants

#endif  // FALCO_CONFIG_HPP_
