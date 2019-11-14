/* falco-diff : compare outputs of two fastqc/falco runs
 *
 * Copyright (C) 2019 Guilherme De Sena Brandine and
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


#ifndef _MODULE_CPP
#define _MODULE_CPP
#include "FastqStats.hpp"
#include "FalcoConfig.hpp"
#include "HtmlMaker.hpp"
#include <string>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <utility>

/* base groups for longer reads, copied from FastQC*/
struct BaseGroup {
  size_t start,end;
  BaseGroup (size_t _start, size_t _end) : start(_start), end(_end) {}
};


class Module {
 private:
  // avoid writing things prior to summarizing
  bool summarized;

 public:
  // The module name displayed in outputs and html
  // GS TODO: automate placing it in html too
  const std::string module_name;

  // the module name lowercased without spaces
  std::string placeholder;

  // the placeholder in the html for the data the module generates
  std::string placeholder_data;

  // the placeholder in the html for the comments
  // (comment start: <!-- and comment end: -->
  std::string placeholder_cs;
  std::string placeholder_ce;

  // placeholder for grade
  std::string placeholder_grade;

  // placeholder for module name
  std::string placeholder_name;

  // pass warn fail
  std::string grade;

  std::string html_data;
  Module(const std::string &_module_name);
  virtual ~Module() = 0;

  /*********************************************/
  /*****Abstract functions to be implemented****/
  /*********************************************/
  // Summarize the module
  virtual void summarize_module(const FastqStats &stats) = 0;

  // Decide if it's a pass/warn/fail
  virtual void make_grade() = 0;

  // write long summary
  virtual void write_module(std::ostream &os) = 0;
  virtual std::string make_html_data() = 0;

  /*********************************************/
  /**************Visible functions**************/
  /*********************************************/
  // Summarizes and registers that it summarized
  void summarize(const FastqStats &stats);

  // Write the module in the FastQC standard, starting with title,
  // pass/warn/fail and then ending with >>END_MODULE
  void write(std::ostream &os);

  // write short summary
  void write_short_summary(std::ostream &os, const std::string &filename);

  // Put html data
  void put_data_on_html(HtmlMaker &html_maker);
};

class ModuleBasicStatistics : public Module {
 private:
  std::string file_type;
  std::string file_encoding;
  std::string filename_stripped;
  size_t avg_read_length;
  size_t avg_gc;
  size_t num_poor;
  size_t min_read_length;
  size_t max_read_length;
  size_t total_sequences;
 public:

  ModuleBasicStatistics(const FalcoConfig &config);
  ~ModuleBasicStatistics() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModulePerBaseSequenceQuality : public Module {
 private:
  // from FastQC: whether to group bases
  bool do_group;
  size_t num_bases;
  size_t num_groups;
  // grade criteria
  size_t base_lower_warn,
         base_lower_error,
         base_median_warn,
         base_median_error;
   size_t num_warn, num_error;
  std::vector<double> group_mean;
  std::vector<size_t> group_ldecile,
                      group_lquartile,
                      group_median,
                      group_uquartile,
                      group_udecile;
  std::vector<BaseGroup> base_groups;

 public:
  ModulePerBaseSequenceQuality(const FalcoConfig &config);
  ~ModulePerBaseSequenceQuality() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModulePerTileSequenceQuality : public Module {
 private:
   double grade_warn, grade_error;
   size_t max_read_length;
   std::unordered_map<size_t, std::vector<double>> tile_position_quality;
   std::vector<size_t> tiles_sorted;
 public:
  ModulePerTileSequenceQuality(const FalcoConfig &config);
  ~ModulePerTileSequenceQuality() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModulePerSequenceQualityScores : public Module {
 private:
   size_t mode_val;
   size_t mode_ind;
   std::array<size_t, FastqStats::kNumQualityValues> quality_count;
   // grade criteria
   size_t mode_warn;
   size_t mode_error;
 public:
  ModulePerSequenceQualityScores(const FalcoConfig &config);
  ~ModulePerSequenceQualityScores() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModulePerBaseSequenceContent : public Module {
 private:
   std::vector<double> a_pct, c_pct, t_pct, g_pct;
   double max_diff;
   size_t num_bases;

   // flag as to whether or not dataset is WGBS
   bool is_bisulfite;
   // for grade
   double sequence_error, sequence_warn;
 public:
  ModulePerBaseSequenceContent(const FalcoConfig &config);
  ~ModulePerBaseSequenceContent() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModulePerSequenceGCContent : public Module {
 private:
   double gc_warn, gc_error;
   double gc_deviation;
   std::array<double, 101> gc_count;
   std::array<double, 101> theoretical_gc_count;

 public:
  ModulePerSequenceGCContent(const FalcoConfig &config);
  ~ModulePerSequenceGCContent() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModulePerBaseNContent : public Module {
 private:
   size_t num_bases;
   // for grade
   size_t grade_n_warn;
   size_t grade_n_error;

   double gc_deviation;
   std::array<size_t, 101> gc_count;
   std::array<size_t, 101> theoretical_gc_count;
   std::vector<double> n_pct;
   // grade vars
   size_t gc_warn, gc_error;
 public:
  ModulePerBaseNContent(const FalcoConfig &config);
  ~ModulePerBaseNContent() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModuleSequenceLengthDistribution : public Module {
 private:
  bool do_grade_error;
  bool do_grade_warn;
  size_t max_read_length;
  std::vector<size_t> sequence_lengths;

  // warn and fail criteria
  bool has_empty_read;
  bool is_all_same_length;
 public:
  ModuleSequenceLengthDistribution(const FalcoConfig &config);
  ~ModuleSequenceLengthDistribution() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModuleSequenceDuplicationLevels : public Module {
 private:
   double seq_total, seq_dedup;

   double grade_dup_warn;
   double grade_dup_error;
   double total_deduplicated_pct;
   std::array<double, 16> percentage_deduplicated;
   std::array<double, 16> percentage_total;
   std::unordered_map<size_t,size_t> counts_by_freq;
 public:
  ModuleSequenceDuplicationLevels(const FalcoConfig &config);
  ~ModuleSequenceDuplicationLevels() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModuleOverrepresentedSequences : public Module {
 private:
  size_t num_reads;
  std::vector<std::pair<std::string, size_t>> overrep_sequences;
  double grade_warn, grade_error;
  const double min_fraction_to_overrepresented = 0.001;
  std::vector<std::pair<std::string, std::string> > contaminants;

  // Function to find the matching contaminant within the list
  std::string get_matching_contaminant(const std::string &seq);
 public:
  ModuleOverrepresentedSequences(const FalcoConfig &config);
  ~ModuleOverrepresentedSequences() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModuleAdapterContent : public Module {
 private:
   // Number of adapters to test
   size_t num_adapters;

   // number of bases to report
   size_t num_bases;

   // Information from config
   std::vector<std::string> adapter_names;
   std::vector<std::string> adapter_seqs;
   std::vector<size_t> adapter_hashes;

   // vector to be reported
   std::vector<std::vector<double>> adapter_pos_pct;
   // minimum percentages for warn/fail
   double grade_warn, grade_error;

   // Aux function to count adapter in a position
   double count_adapter (const std::vector<size_t> &kmer_count,
                         const size_t pos,
                         const size_t adapter_hash,
                         const size_t adapter_size,
                         const size_t kmer_size);
 public:
  ModuleAdapterContent(const FalcoConfig &config);
  ~ModuleAdapterContent() {}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};

class ModuleKmerContent : public Module {
 private:
   size_t num_kmer_bases;
   size_t kmer_size;
   size_t num_kmers;
   size_t num_seen_kmers;

   double grade_warn, grade_error;
   std::array<size_t, FastqStats::kNumBases> pos_kmer_count;
   std::vector<size_t> total_kmer_counts;
   std::vector<double> obs_exp_max;
   std::vector<size_t> where_obs_exp_is_max;
   std::vector<std::pair<size_t, double>> kmers_to_report;
 public:
  ModuleKmerContent(const FalcoConfig &config);
  ~ModuleKmerContent(){}
  void summarize_module(const FastqStats &stats);
  void make_grade();
  void write_module(std::ostream &os);
  std::string make_html_data();
};
#endif

