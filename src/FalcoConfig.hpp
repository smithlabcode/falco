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

#ifndef FALCO_CONFIG_HPP
#define FALCO_CONFIG_HPP

#include <string>
#include <vector>
#include <unordered_map>
#include <utility>

#include "aux.hpp"


/*************************************************************
 ******************** CUSTOM CONFIGURATION *******************
 *************************************************************/

// config from options, constants, magic numbers, etc
struct FalcoConfig {
  static const std::string FalcoVersion;
  FalcoConfig();  // set magic defaults

  /************************************************************
   *************** FASTQC OPTION PARSER************************
   ************************************************************/
  bool casava;  // files from raw casava output
  bool nanopore;  // fast5 format
  bool nofilter;  // if running with --casava flag
  bool extract;  // if set the zipped file will be uncompressed
  bool nogroup;  // disable grouping of bases for reads >50bp
  bool compressed;  // whether or not to inflate file
  bool quiet;
  size_t min_length;  // lower limit in sequence length to be shown in report
  size_t threads;  // number of threads to read multiple files in parallel
  std::string format;  // force file format
  std::string contaminants_file;  // custom contaminants file
  std::string adapters_file;  // adapters file
  std::string limits_file;  // file with limits and options and custom analyses
  static const std::string html_template; // the html for the template
  std::string tmpdir;  // dir for temp files when generating report images

  // config on how to handle reads
  bool do_duplication,
       do_kmer,
       do_n_content,
       do_overrepresented,
       do_quality_base,
       do_sequence,
       do_gc_sequence,
       do_quality_sequence,
       do_tile,
       do_adapter,
       do_adapter_optimized,
       do_sequence_length;

  /************************************************************
   *************** FASTQC LIMITS *******************************
   ************************************************************/
  // These will become const bools in the stream reader
  std::unordered_map<std::string,
                     std::unordered_map<std::string, double> > limits;
  static const std::vector<std::string> values_to_check;

  /*************** CONTAMINANTS *****************/
  // below: first = name, scond = seq
  std::vector<std::pair<std::string, std::string> > contaminants;

  /*************** ADAPTERS *********************/
  // Name (eg: Illumina Small RNA adapter)
  std::vector<std::string> adapter_names;

  // Actual string sequence (eg: ATTGCCACA)
  std::vector<std::string> adapter_seqs;

  // two-bit hash of the sequence above
  std::vector<size_t> adapter_hashes;

  size_t adapter_size;
  size_t shortest_adapter_size;
  /************************************************************
   ******* ADDITIONAL INFORMATION ABOUT THE SAMPLE ************
   ************************************************************/
  bool is_bisulfite;
  bool is_reverse_complement;
 
  /*************** DEFINE FILE TYPE ************/

  // IO
  bool is_sam, is_bam, is_fastq, is_fastq_gz;
  std::string filename;
  std::string filename_stripped;

  /*********** FUNCTIONS TO READ FILES *************/
  void define_file_format();
  void read_limits();  // populate limits hash map
  void read_adapters();
  void read_contaminants_file();

  void setup();
};

/*************************************************************
 ******************** ALL MAGIC NUMBERS **********************
 *************************************************************/
namespace Constants {
  // log of a power of two, to use in bit shifting for fast index acces
  // returns the log2 of a number if it is a power of two, or zero
  // otherwise
  constexpr size_t
  log2exact(size_t v) {
    return (63 -
            ((v & 0x00000000FFFFFFFF) ? 32 : 0) -
            ((v & 0x0000FFFF0000FFFF) ? 16 : 0) -
            ((v & 0x00FF00FF00FF00FF) ?  8 : 0) -
            ((v & 0x0F0F0F0F0F0F0F0F) ?  4 : 0) -
            ((v & 0x3333333333333333) ?  2 : 0) -
            ((v & 0x5555555555555555) ?  1 : 0));
  }

  const size_t kmer_size = 7;
  const size_t max_adapters = 128;

  // number of bases for static allocation.
  const size_t num_static_bases = 500;

  // Value to subtract quality characters to get the actual quality value
  const size_t quality_zero = 33;  // The ascii for the lowest quality

  // Smallest power of two that comprises all possible Illumina quality values.
  // Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of
  // is to avoid double pointer jumps and to get indices with bit shifts.
  const size_t num_quality_values = 128;

  // How many possible nucleotides (must be power of 2!)
  const size_t num_nucleotides = 4;  // A = 00,C = 01,T = 10,G = 11

  /************* DUPLICATION ESTIMATES *************/
  // Number of unique sequences to see before stopping counting sequences
  const size_t unique_reads_stop_counting = 1e5;

  // Maximum read length to store the entire read in memory
  const size_t unique_reads_max_length = 75;

  // Prefix size to cut if read length exceeds the value above
  const size_t unique_reads_truncate = 50;

  /****Bit shifts as instructions for the std::arrays***/
  // for matrices that count stats per nucleotide
  const size_t bit_shift_base = log2exact(num_nucleotides);

  // for matrices that count stats for quality value
  const size_t bit_shift_quality = log2exact(num_quality_values);

  // bit shift for adapters, log(128) = 7
  const size_t bit_shift_adapter = log2exact(max_adapters);

  // we shift 14 bits when reading a kmer, two bits per base
  const size_t bit_shift_kmer = bit_shift_base*kmer_size;

  // mask to get only the first 2*k bits of the sliding window
  const size_t kmer_mask = (1ull << (bit_shift_kmer)) - 1;
};

#endif
