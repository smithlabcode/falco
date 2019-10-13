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

#ifndef FASTQSTATS_HPP
#define FASTQSTATS_HPP

#include <string>
#include <vector>
#include <iostream>
#include <unordered_map>
#include "FalcoConfig.hpp"

// log of qa power of two, to use in bit shifting for fast index acces
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

/*************************************************************
 ******************** FASTQ STATS ****************************
 *************************************************************/

struct FastqStats {
  // number of bases for static allocation.
  static const size_t kNumBases = 500;

  // Value to subtract quality characters to get the actual quality value
  static const size_t kBaseQuality = 33;  // The ascii for the lowest quality

  // Smallest power of two that comprises all possible Illumina quality values.
  // Illumina gives qualities from 0 to 40, therefore we set it as 64. Power of
  // is to avoid double pointer jumps and to get indices with bit shifts.
  static const size_t kNumQualityValues = 64;

  // How many possible nucleotides (must be power of 2!)
  static const size_t kNumNucleotides = 4;  // A = 00,C = 01,T = 10,G = 11

  // maximum tile value
  static const size_t kMinTileNumber = 10000;

  // Maximum number of bases for which to do kmer statistics
  static const size_t kKmerMaxBases = 500;

  /************* DUPLICATION ESTIMATES *************/
  // Number of unique sequences to see before stopping counting sequences
  static const size_t kDupUniqueCutoff = 1e5;

  // Maximum read length to store the entire read in memory
  static const size_t kDupReadMaxSize = 75;

  // Prefix size to cut if read length exceeds the value above
  static const size_t kDupReadTruncateSize = 50;

  // Kmer size given as input
  static const size_t kmer_size = 7;

  // Bit shifts as instructions for the std::arrays
  static const size_t kBitShiftNucleotide = log2exact(kNumNucleotides);
  static const size_t kBitShiftQuality = log2exact(kNumQualityValues);
  static const size_t kBitShiftKmer = 2 * kmer_size;  // two bits per base

 public:
  /*********** SINGLE NUMBERS FROM THE ENTIRE FASTQ ****************/
  // Number of unique sequences seen thus far
  size_t num_unique_seen;

  // How many reads were processed before num_unique_seen = kDupUniqueCutoff
  size_t count_at_limit;

  size_t total_bases;  // sum of all bases in all reads
  size_t avg_read_length;  // average of all read lengths
  size_t num_reads;  // total number of lines read
  size_t num_reads_kmer;  // number of reads in which kmer counts was performed
  size_t min_read_length;  // minimum read length seen
  size_t max_read_length;  // total number of lines read
  size_t num_poor;  // reads whose average quality was <= poor
  size_t num_extra_bases;  // number of bases outside of buffer

  // mask to get only the first 2*k bits of the sliding window
  size_t kmer_mask;
  size_t total_gc;
  double avg_gc;  // (sum of g bases + c bases) / (num_reads)
  double total_deduplicated_pct;  // number of reads left if deduplicated

  /*********************************************************
   *********** METRICS COLLECTED DURING IO *****************
   *********************************************************/
  /*********** PER BASE METRICS ****************/

  // counts the number of bases in every read position
  std::array<size_t, kNumNucleotides * kNumBases> base_count;  // ATGC
  std::array<size_t, kNumNucleotides * kNumBases> n_base_count;  // N

  /*********** PER QUALITY VALUE METRICS ****************/
  // Counts of quality in each base position
  std::array<size_t, kNumQualityValues * kNumBases> position_quality_count;

  // Counts of average quality (truncated) per sequence
  std::array<size_t, kNumQualityValues> quality_count;

  /*********** PER GC VALUE METRICS ****************/
  // histogram of GC fraction in each read from 0 to 100%
  std::array<size_t, 101> gc_count;
  std::array<double, 101> smooth_gc_count;
  std::array<double, 101> theoretical_gc_count;

    /*********** PER READ METRICS ***************/
  // Distribution of read lengths
  std::array<size_t, kNumBases> read_length_freq;
  std::array<size_t, kNumBases> cumulative_read_length_freq;

  /*********** PER TILE SEQUENCE QUALITY OVERSERQUENCES ********/
  std::unordered_map <size_t, std::vector<double> > tile_position_quality;
  std::unordered_map <size_t, std::vector<size_t> > tile_position_count;
  std::unordered_map <size_t, size_t> tile_count;

  /*********** SUMMARY ARRAYS **************/
  // Quantiles for the position_quality_count
  std::array<size_t, kNumBases> ldecile, lquartile, median, uquartile, udecile;
  std::array<double, kNumBases> mean;

  // For sequence duplication levels
  // 1 to 9, >10, >50, >100, >500, >1k, >5k, >10k+
  std::array<double, 16> percentage_deduplicated;
  std::array<double, 16> percentage_total;

  // Percentages for per base sequence content
  std::array<double, kNumBases> a_pct,
                           c_pct,
                           t_pct,
                           g_pct,
                           n_pct;

  /*********** SLOW STUFF *******************/
  // Leftover memory using dynamic allocation
  std::vector<size_t> long_base_count;
  std::vector<size_t> long_n_base_count;
  std::vector<size_t> long_position_quality_count;
  std::vector<size_t> long_read_length_freq;
  std::vector<size_t> long_cumulative_read_length_freq;
  std::vector<size_t> long_ldecile, long_lquartile, long_median,
                 long_uquartile, long_udecile;
  std::vector<double> long_mean;
  std::vector<double> long_a_pct,
                 long_c_pct,
                 long_t_pct,
                 long_g_pct,
                 long_n_pct;


  /********** KMER FREQUENCY ****************/
  // A 2^K + 1 std::vector to count all possible kmers
  std::vector<size_t> kmer_count;

  /********** ADAPTER COUNT *****************/
  // Percentage of times we saw each adapter in each position
  std::unordered_map <size_t, std::vector <double>> kmer_by_base;
  /*********** DUPLICATION ******************/
  std::unordered_map <std::string, size_t> sequence_count;

  /*********** OVERREPRESENTED SERQUENCES ********/
  std::vector <std::pair<std::string, size_t>> overrep_sequences;

  /*********************************************************
   *********** METRICS SUMMARIZED AFTER IO *****************
   *********************************************************/

  // I need this to know what to divide each base by
  // when averaging content, bases, etc. It stores, for every element i, how
  // many reads are of length >= i, ie, how many reads actually have a
  // nucleotide at position i

  /*********** PASS WARN FAIL MESSAGE FOR EACH METRIC **************/
  std::string pass_basic_statistics,
         pass_per_base_sequence_quality,
         pass_per_tile_sequence_quality,
         pass_per_sequence_quality_scores,
         pass_per_base_sequence_content,
         pass_per_sequence_gc_content,
         pass_per_base_n_content,
         pass_sequence_length_distribution,
         pass_overrepresented_sequences,
         pass_duplicate_sequences,
         pass_adapter_content;

    /**************** FUNCTIONS ****************************/

  // Default constructor that zeros everything
  FastqStats();

  // Allocation of more read positions
  void allocate_new_base(const bool ignore_tile);

  /******* DUPLICATION AND OVERREPRESENTATION *******/
  // Makes a hash map with keys as 32-bit suffixes and values as all the
  // candidate frequent sequences with that given suffix

  // Summarize all statistics we need before writing
  void summarize(FalcoConfig &config);

  // Writes equivalent of fastqc_data.txt
  void write(std::ostream &os, const FalcoConfig &config);

  // Writes equivalent of fastqc's summary.txt
  void write_summary(std::ostream &os, const FalcoConfig &config);
};


#endif
