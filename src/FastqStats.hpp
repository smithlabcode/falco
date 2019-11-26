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
#include <array>
#include "FalcoConfig.hpp"
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

/******************* BEGIN COPY FROM FASTQC *****************/
struct GCModelValue {
  int percent;
  double increment;
  GCModelValue() {
    percent = 0;
    increment = 0.0;
  }
};
struct GCModel {
  std::vector<std::vector<GCModelValue>> models;
  GCModel() {}
  GCModel(const int read_length) {
    if (read_length == 0)
      return;

    // Number of counts that goes into each bin
    std::array<size_t, 101> claiming_counts;
    claiming_counts.fill(0);


    // Iterate over all possible gc counts (pos)
    for (int pos = 0; pos <= read_length; pos++) {
      double low_count = static_cast<double>(pos) - 0.5;
      double high_count = static_cast<double>(pos) + 0.5;

      if (low_count < 0) low_count = 0;
      if (high_count < 0) high_count = 0;
      if (high_count > read_length) high_count = read_length;
      if (low_count > read_length) low_count = read_length;

      int low_pct = (int)round(100 * low_count / 
                    static_cast<double>(read_length));
      int high_pct = (int)round(100 * high_count / 
                    static_cast<double>(read_length));

      for(int p = low_pct; p <= high_pct; p++) {
        claiming_counts[p]++;
      }
    }

    // We now do a second pass to make up the model using the weightings
    // we calculated previously.
    for (int pos = 0; pos <= read_length; pos++) {
      double low_count = static_cast<double>(pos) - 0.5;
      double high_count = static_cast<double>(pos) + 0.5;

      if (low_count < 0) low_count = 0;
      if (high_count < 0) high_count = 0;
      if (high_count > read_length) high_count = read_length;
      if (low_count > read_length) low_count = read_length;

      // Check the bins in which percentages must be put
      int low_pct = (int)round((100 * low_count) /
                      static_cast<double>(read_length));

      int high_pct = (int)round((100 * high_count) /
                      static_cast<double>(read_length));

      // Add a new vector of values
      models.push_back(
        std::vector<GCModelValue>(high_pct - low_pct + 1, GCModelValue())
      );

      // populates the increment in each bin
      for (int p = low_pct; p <= high_pct; ++p) {
        models[pos][p - low_pct].percent = p;
        models[pos][p - low_pct].increment =
          1.0 / static_cast<double>(claiming_counts[p]);
      }
    }
  }
};

/********************** END COPY FROM FASTQC *************/

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
  static const size_t kNumQualityValues = 128;

  // How many possible nucleotides (must be power of 2!)
  static const size_t kNumNucleotides = 4;  // A = 00,C = 01,T = 10,G = 11

  /************* DUPLICATION ESTIMATES *************/
  // Number of unique sequences to see before stopping counting sequences
  static const size_t kDupUniqueCutoff = 1e5;

  // Maximum read length to store the entire read in memory
  static const size_t kDupReadMaxSize = 75;

  // Prefix size to cut if read length exceeds the value above
  static const size_t kDupReadTruncateSize = 50;

   // Bit shifts as instructions for the std::arrays
  static const size_t kBitShiftNucleotide = log2exact(kNumNucleotides);
  static const size_t kBitShiftQuality = log2exact(kNumQualityValues);

  /************ KMER CONSTANTS **********/
  // Kmer size given as input
  static const size_t kmer_size = 7;

  // we shift 14 bits when reading a kmer, two bits per base
  static const size_t kBitShiftKmer = 2 * kmer_size;

  // mask to get only the first 2*k bits of the sliding window
  static const size_t kmer_mask = (1ll << (2*kmer_size)) - 1;

  /************ ADAPTER CONSTANTS **********/
  // Maximum number of adapters
  static const size_t max_adapters = 128;

  // bit shift for adapters, log(100) = 7
  static const size_t kBitShiftAdapter = log2exact(max_adapters);


 public:
  /*********** SINGLE NUMBERS FROM THE ENTIRE FASTQ ****************/
  // Number of unique sequences seen thus far
  size_t num_unique_seen;

  // How many reads were processed before num_unique_seen = kDupUniqueCutoff
  size_t count_at_limit;

  size_t total_bases;  // sum of all bases in all reads
  size_t num_reads;  // total number of lines read
  size_t min_read_length;  // minimum read length seen
  size_t max_read_length;  // total number of lines read
  size_t num_poor;  // reads whose average quality was <= poor
  size_t num_extra_bases;  // number of bases outside of buffer
  size_t total_gc; // sum of all G+C bases in all reads

  // Pre-calculated GC model increments
  static const std::array<GCModel, kNumBases> gc_models;

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
  std::array<double, 101> gc_count;

    /*********** PER READ METRICS ***************/
  // Distribution of read lengths
  std::array<size_t, kNumBases> read_length_freq;
  std::array<size_t, kNumBases> cumulative_read_length_freq;

  /*********** PER TILE SEQUENCE QUALITY OVERSERQUENCES ********/
  std::unordered_map <size_t, std::vector<double> > tile_position_quality;
  std::unordered_map <size_t, std::vector<size_t> > tile_position_count;

  /*********** SLOW DATA STRUCTURES FOR LONGER READS ************/
  // Leftover memory using dynamic allocation
  std::vector<size_t> long_base_count;
  std::vector<size_t> long_n_base_count;
  std::vector<size_t> long_position_quality_count;
  std::vector<size_t> long_read_length_freq;
  std::vector<size_t> long_cumulative_read_length_freq;

  /********** KMER FREQUENCY ****************/
  // A (4^K + 1)*kNumBases std::vector to count all possible kmers
  std::vector<size_t> kmer_count;

  // How many kmers were counted in each position
  std::array<size_t, kNumBases> pos_kmer_count;

  // How many adapters were counted in each position
  std::array<size_t, kNumBases> pos_adapter_count;

  /*********** DUPLICATION ******************/
  // First 100k unique sequences and how often they were seen
  std::unordered_map <std::string, size_t> sequence_count;

  /**************** FUNCTIONS ****************************/
  // Default constructor that zeros everything
  FastqStats();

  // Allocation of more read positions
  void allocate_new_base(const bool ignore_tile);

  void summarize();

  // Given an input fastqc_data.txt file, populate the statistics with it
  void read(std::istream &is);
};
#endif
