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

#include "FastqStats.hpp"

#include <algorithm>
#include <iomanip>
#include <cmath>
using std::string;
using std::vector;
using std::array;
using std::unordered_map;
using std::sort;
using std::min;
using std::max;
using std::ostream;
using std::pair;
using std::transform;
using std::toupper;
using std::setprecision;

// To make the gc models static const
static array<GCModel, FastqStats::kNumBases>
make_gc_models () {
  array<GCModel, FastqStats::kNumBases> ans;
  for (size_t i = 0; i < FastqStats::kNumBases; ++i) {
    ans[i] = GCModel(i);
  }
  return ans;
}

/****************************************************************/
/******************** FASTQ STATS FUNCTIONS *********************/
/****************************************************************/
// Default constructor
FastqStats::FastqStats() {
  lowest_char = std::numeric_limits<char>::max();
  encoding_offset = 0;

  total_bases = 0;
  num_extra_bases = 0;
  total_gc = 0;
  num_reads = 0;
  empty_reads = 0;
  min_read_length = 0;
  max_read_length = 0;
  num_poor = 0;

  num_unique_seen = 0;
  count_at_limit = 0;

  // Initialize IO arrays
  base_count.fill(0);
  n_base_count.fill(0);
  read_length_freq.fill(0);
  quality_count.fill(0);
  gc_count.fill(0);
  position_quality_count.fill(0);
  pos_kmer_count.fill(0);
  pos_adapter_count.fill(0);
  kmer_count = vector<size_t>(kNumBases*(Constants::kmer_mask + 1), 0);
}

// Initialize as many gc models as fast bases
const array<GCModel, FastqStats::kNumBases>
FastqStats::gc_models = make_gc_models();

// When we read new bases, dynamically allocate new space for their statistics
void
FastqStats::allocate_new_base(const bool ignore_tile) {
  for (size_t i = 0; i < kNumNucleotides; ++i) {
    long_base_count.push_back(0);
  }

  long_n_base_count.push_back(0);

  // space for quality boxplot
  for (size_t i = 0; i < kNumQualityValues; ++i)
    long_position_quality_count.push_back(0);

  long_read_length_freq.push_back(0);

  // space for tile quality in each position.
  //
  if (!ignore_tile) {
    for (auto  &v : tile_position_quality) {
      v.second.push_back(0);
    }
    for (auto  &v : tile_position_count) {
      v.second.push_back(0);
    }
  }

  // Successfully allocated space for a new base
  ++num_extra_bases;
}

// Calculates all summary statistics and pass warn fails
void
FastqStats::summarize() {
  // Cumulative read length frequency
  size_t cumulative_sum = 0;
  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases) {
      cumulative_sum += read_length_freq[i];
      if (read_length_freq[i] > 0)
        if (min_read_length == 0)
          min_read_length = i + 1;
    }
    else {
      cumulative_sum += long_read_length_freq[i - kNumBases];
      if (long_read_length_freq[i - kNumBases] > 0)
        if (min_read_length == 0)
          min_read_length = i + 1;
    }
  }

  for (size_t i = 0; i < max_read_length; ++i) {
    if (i < kNumBases) {
      cumulative_read_length_freq[i] = cumulative_sum;
      cumulative_sum -= read_length_freq[i];
    }
    else {
      long_cumulative_read_length_freq.push_back(cumulative_sum);
      cumulative_sum -= long_read_length_freq[i - kNumBases];
    }
  }
}

