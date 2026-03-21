/*
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
#include "FalcoConfig.hpp"
#include "FastqStats.hpp"
#include "HtmlMaker.hpp"
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

/* base groups for longer reads, copied from FastQC*/
struct BaseGroup {
  std::size_t start{};
  std::size_t end{};
};

struct Module {
  const std::string module_name{};

  // avoid writing things prior to summarizing
  bool summarized{};

  // The module name displayed in outputs and html
  // GS TODO: automate placing it in html too

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
  Module(const std::string &module_name);
  virtual ~Module() = 0;

  // functions to be in child classes

  // Summarize the module
  virtual void
  summarize_module(FastqStats &stats) = 0;

  // Decide if it's a pass/warn/fail
  virtual void
  make_grade() = 0;

  // write long summary
  virtual void
  write_module(std::ostream &os) = 0;
  virtual std::string
  make_html_data() = 0;

  // Summarizes and registers that it summarized
  void
  summarize(FastqStats &stats);

  // Write the module in the FastQC standard, starting with title,
  // pass/warn/fail and then ending with >>END_MODULE
  void
  write(std::ostream &os);

  // write short summary
  void
  write_short_summary(std::ostream &os, const std::string &filename);

  // Put html data
  void
  put_data_on_html(HtmlMaker &html_maker);
};

struct ModuleBasicStatistics : public Module {
  static const std::string module_name;
  bool is_nanopore{};
  std::string file_type;
  std::string file_encoding;
  std::string filename_stripped;
  std::size_t avg_read_length{};
  double avg_gc{};
  std::size_t num_poor{};
  std::size_t min_read_length{};
  std::size_t max_read_length{};
  std::size_t total_bases{};
  std::size_t total_sequences{};
  ModuleBasicStatistics(const FalcoConfig &config);
  ~ModuleBasicStatistics() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
  void
  read_data_line(const std::string &line);
};

struct ModulePerBaseSequenceQuality : public Module {
  static const std::string module_name;
  // from FastQC: whether to group bases
  bool do_group{};
  std::size_t num_bases{};
  std::size_t num_groups{};
  std::size_t base_lower_warn{};   // grade criteria
  std::size_t base_lower_error{};  // grade criteria
  std::size_t base_median_warn{};
  std::size_t base_median_error{};
  std::size_t num_warn{};
  std::size_t num_error{};
  std::vector<double> group_mean;
  std::vector<double> group_ldecile;
  std::vector<double> group_lquartile;
  std::vector<double> group_median;
  std::vector<double> group_uquartile;
  std::vector<double> group_udecile;
  std::vector<BaseGroup> base_groups;

  ModulePerBaseSequenceQuality(const FalcoConfig &config);
  ~ModulePerBaseSequenceQuality() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  void
  read_data_line([[maybe_unused]] const std::string &line);
  std::string
  make_html_data();
};

struct ModulePerTileSequenceQuality : public Module {
  double grade_warn{};
  double grade_error{};
  std::size_t max_read_length{};
  std::unordered_map<std::size_t, std::vector<double>> tile_position_quality;
  std::vector<std::size_t> tiles_sorted;

  static const std::string module_name;
  ModulePerTileSequenceQuality(const FalcoConfig &config);
  ~ModulePerTileSequenceQuality() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModulePerSequenceQualityScores : public Module {
  std::size_t mode_val{};
  std::size_t mode_ind{};
  std::size_t offset{};
  std::array<std::size_t, FastqStats::kNumQualityValues> quality_count{};
  // grade criteria
  std::size_t mode_warn{};
  std::size_t mode_error{};

  static const std::string module_name;
  ModulePerSequenceQualityScores(const FalcoConfig &config);
  ~ModulePerSequenceQualityScores() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModulePerBaseSequenceContent : public Module {
  static const std::string module_name;
  bool do_group{};
  std::vector<double> a_pct;
  std::vector<double> c_pct;
  std::vector<double> t_pct;
  std::vector<double> g_pct;
  double max_diff{};
  std::size_t num_bases{};

  bool is_bisulfite{};           // flag for dataset is WGBS
  bool is_reverse_complement{};  // indicates test T vs C instead of A vs G

  double sequence_error{};  // for grade
  double sequence_warn{};   // for grade

  std::size_t num_groups{};
  std::vector<BaseGroup> base_groups;

  ModulePerBaseSequenceContent(const FalcoConfig &config);
  ~ModulePerBaseSequenceContent() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModulePerSequenceGCContent : public Module {
  static const std::string module_name;
  double gc_warn{};
  double gc_error{};
  double gc_deviation{};
  std::array<double, 101> gc_count;
  std::array<double, 101> theoretical_gc_count;

  ModulePerSequenceGCContent(const FalcoConfig &config);
  ~ModulePerSequenceGCContent() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModulePerBaseNContent : public Module {
  static const std::string module_name;
  std::size_t num_bases{};
  std::size_t grade_n_warn{};   // for grade
  std::size_t grade_n_error{};  // for grade
  double max_n_pct{};
  std::array<std::size_t, 101> gc_count;
  std::array<std::size_t, 101> theoretical_gc_count;
  std::vector<double> n_pct;
  std::size_t gc_warn{};   // grade vars
  std::size_t gc_error{};  // grade vars
  bool do_group{};
  std::size_t num_groups{};
  std::vector<BaseGroup> base_groups;

  ModulePerBaseNContent(const FalcoConfig &config);
  ~ModulePerBaseNContent() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModuleSequenceLengthDistribution : public Module {
  static const std::string module_name;
  bool do_grade_error{};
  bool do_grade_warn{};
  std::size_t max_read_length{};
  std::vector<std::size_t> sequence_lengths;

  // warn and fail criteria
  bool is_all_same_length{};
  std::size_t empty_reads{};

  bool do_group{};
  std::size_t num_groups{};
  std::vector<BaseGroup> base_groups;

  ModuleSequenceLengthDistribution(const FalcoConfig &config);
  ~ModuleSequenceLengthDistribution() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModuleSequenceDuplicationLevels : public Module {
  static const std::string module_name;
  double seq_total{};
  double seq_dedup{};

  double grade_dup_warn{};
  double grade_dup_error{};
  double total_deduplicated_pct{};
  std::array<double, 16> percentage_deduplicated{};
  std::array<double, 16> percentage_total{};
  std::unordered_map<std::size_t, std::size_t> counts_by_freq;

  ModuleSequenceDuplicationLevels(const FalcoConfig &config);
  ~ModuleSequenceDuplicationLevels() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModuleOverrepresentedSequences : public Module {
  static constexpr auto min_fraction_to_overrepresented = 0.001;
  static const std::string module_name;

  std::size_t num_reads{};
  std::vector<std::pair<std::string, std::size_t>> overrep_sequences;
  double grade_warn{};
  double grade_error{};
  std::vector<std::pair<std::string, std::string>> contaminants;

  // Function to find the matching contaminant within the list
  std::string
  get_matching_contaminant(const std::string &seq);

  ModuleOverrepresentedSequences(const FalcoConfig &config);
  ~ModuleOverrepresentedSequences() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModuleAdapterContent : public Module {
  static const std::string module_name;
  std::size_t num_adapters{};  // Number of adapters to test
  std::size_t num_bases{};     // number of bases to report
  // adapter size to know how many bases to report
  std::size_t adapter_size{};

  // info from config
  std::vector<std::string> adapter_names;
  std::vector<std::string> adapter_seqs;
  std::vector<std::size_t> adapter_hashes;
  std::size_t shortest_adapter_size{};

  // to be reported
  std::vector<std::vector<double>> adapter_pos_pct;
  // min minimum cutoffs for warn/fail (percentages)
  double grade_warn{};
  double grade_error{};

  // Aux function to count adapter in a position
  double
  count_adapter(const std::vector<std::size_t> &kmer_count,
                const std::size_t pos, const std::size_t adapter_hash,
                const std::size_t adapter_size, const std::size_t kmer_size);

  ModuleAdapterContent(const FalcoConfig &config);
  ~ModuleAdapterContent() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};

struct ModuleKmerContent : public Module {
  static constexpr std::size_t MIN_OBS_EXP_TO_REPORT = 5;
  static constexpr std::size_t MAX_KMERS_TO_REPORT = 20;
  static constexpr std::size_t MAX_KMERS_TO_PLOT = 10;
  static const std::string module_name;

  std::size_t num_kmer_bases{};
  std::size_t kmer_size{};
  std::size_t num_kmers{};
  std::size_t num_seen_kmers{};

  double grade_warn{};
  double grade_error{};
  std::array<std::size_t, FastqStats::SHORT_READ_THRESHOLD> pos_kmer_count{};
  std::vector<std::size_t> total_kmer_counts;
  std::vector<double> obs_exp_max;
  std::vector<std::size_t> where_obs_exp_is_max;
  std::vector<std::pair<std::size_t, double>> kmers_to_report;

  ModuleKmerContent(const FalcoConfig &config);
  ~ModuleKmerContent() {}
  void
  summarize_module(FastqStats &stats);
  void
  make_grade();
  void
  write_module(std::ostream &os);
  std::string
  make_html_data();
};
#endif
