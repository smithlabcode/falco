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

  FalcoConfig();  // set magic defaults

  /************************************************************
   *************** MY UNIVERSAL CONSTANTS *********************
   ************************************************************/
  // threshold for a sequence to be considered  poor quality
  size_t kPoorQualityThreshold;

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
  size_t kmer_size;  // kmer size
  std::string format;  // force file format
  std::string contaminants_file;  // custom contaminants file
  std::string adapters_file;  // adapters file
  std::string limits_file;  // file with limits and options and custom analyses
  std::string html_file;  // file with limits and options and custom analyses
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
#endif
