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

#include "FalcoConfig.hpp"

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

using std::string;
using std::vector;
using std::unordered_map;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::runtime_error;
using std::istringstream;

/********************************************************************/
/**************************** AUX FUNCTIONS *************************/
/********************************************************************/

// Check if line is not a comment or newline
inline bool
is_content_line (const string &line) {
  // comment
  if (line[0] == '#')
    return false;

  // newline
  if (line.size() <= 1)
    return false;

  return true;
}

// Check existance of config files
inline bool
file_exists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

// Check if a std::string ends with another,
// to be use to figure out the file format
inline bool
endswith(std::string const &value, std::string const &ending) {
  if (ending.size() > value.size()) {
    return false;
  }
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

// Removes absolute path from a file
static string
strip_path(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  return full_path.substr(start);
}

/********************************************************************/
/******************** FALCOCONFIG FUNCTIONS *************************/
/********************************************************************/
// Sets magic numbers
FalcoConfig::FalcoConfig() {
  kPoorQualityThreshold = 20;
  casava = false;
  nanopore = false;
  nofilter = false;
  extract = false;
  nogroup = false;
  is_bisulfite = false;
  min_length = 0;
  format = "";
  threads = 1;
  contaminants_file = string(PROGRAM_PATH) +
                      "/Configuration/contaminant_list.txt";
  adapters_file = string(PROGRAM_PATH) +
                      "/Configuration/adapter_list.txt";
  limits_file = string(PROGRAM_PATH) +
                      "/Configuration/limits.txt";
  html_file = string(PROGRAM_PATH) +
                      "/Configuration/template.html";

  // Template has to be checked here
  if (!file_exists(html_file)) {
    throw runtime_error("HTML template file does not exist: " + html_file);
  }
  kmer_size = 7;
  quiet = false;
  tmpdir = ".";

  is_sam = false;
  is_bam = false;
  is_fastq = false;
  is_fastq_gz = false;
}

const vector<string> FalcoConfig::values_to_check({
    "duplication",
    "kmer",
    "n_content",
    "overrepresented",
    "quality_base",
    "sequence",
    "gc_sequence",
    "quality_sequence",
    "tile",
    "sequence_length",
    "adapter",
    "duplication",
    "kmer",
    "n_content",
    "overrepresented",
    "quality_base_lower",
    "quality_base_median",
    "sequence",
    "gc_sequence",
    "quality_sequence",
    "tile",
    "sequence_length",
    "adapter"
  });

void
FalcoConfig::setup() {
  // Now check for the file format (FASTQ/SAM/BAM, compressed or not)
  define_file_format();

  // Get filename without absolute path
  filename_stripped = strip_path(filename);

  // read which modules to run and the cutoffs for pass/warn/fail
  read_limits();

  // Read files for appropriate modules
  if (limits["adapter"]["ignore"] == 0.0)
    read_adapters();
  if (limits["adapter"]["ignore"] == 0.0)
    read_contaminants_file();
}

void
FalcoConfig::define_file_format() {
  if (format == "") {
    if (endswith(filename, "sam")) {
      is_sam = true;
    }
    if (endswith(filename, "bam")) {
      is_bam = true;
    }
    if (endswith(filename, "fastq.gz")) {
      is_fastq_gz = true;
    }
    if (endswith(filename, "fq.gz")) {
      is_fastq_gz = true;
    }
    if (endswith(filename, "fastq")) {
      is_fastq = true;
    }
    if (endswith(filename, "fq")) {
      is_fastq = true;
    }
  }
}

void
FalcoConfig::read_limits() {

  ifstream in(limits_file);
  if (!in)
    throw runtime_error("limits file does not exist: " + limits_file);

  // Variables to parse lines
  string line, instruction;
  double value;
  while (getline(in, line)) {
    // Checks if the line has something to be parsed
    if (is_content_line (line)) {
      istringstream iss(line);

      // Every line is a limit, warn/error/ignore and the value
      string limit;
      if (!(iss >> limit >> instruction >> value))
        throw runtime_error("malformed limits line: \"" + line + "\"");

      if (find(begin(values_to_check), end(values_to_check), limit)
          == end(values_to_check))
        throw runtime_error("unknown limit option: " + limit);

      if (instruction != "warn" &&
          instruction != "error" &&
          instruction != "ignore")
        throw runtime_error("unknown instruction for limit " +
                            limit + ": " + instruction);

      limits[limit][instruction] = value;
    }
  }

  for (auto v : values_to_check)
    if (limits.count(v) == 0)
      throw runtime_error("instruction for limit " + v +
                          " not found in file " + limits_file);

  // Get useful data from config that tells us which analyses to skip
  do_duplication = (limits["duplication"]["ignore"] == 0.0);
  do_kmer = (limits["kmer"]["ignore"] == 0.0);
  do_n_content = (limits["n_content"]["ignore"] == 0.0);
  do_overrepresented = (limits["overrepresented"]["ignore"] == 0.0);
  do_quality_base = (limits["quality_base"]["ignore"] == 0.0);
  do_sequence = (limits["sequence"]["ignore"] == 0.0);
  do_gc_sequence = (limits["gc_sequence"]["ignore"] == 0.0);
  do_quality_sequence= (limits["quality_sequence"]["ignore"] == 0.0);
  do_tile = (limits["tile"]["ignore"] == 0.0);
  do_adapter = (limits["adapter"]["ignore"] == 0.0);
  do_sequence_length = (limits["sequence_length"]["ignore"] == 0.0);
}

void
FalcoConfig::read_adapters() {
  ifstream in(adapters_file);
  if (!in)
    throw runtime_error("adapter file not found: " + adapters_file);

  string line, _tmp;
  vector<string> line_by_space;
  string adapter_name, adapter_seq;
  size_t adapter_hash;

  // The adapters file has a space separated name, and the last instance is
  // the biological sequence
  while (getline(in, line)) {
    if (is_content_line(line)) {
      adapter_name = "";
      adapter_seq = "";
      istringstream iss(line);
      while (iss >> _tmp) {
        line_by_space.push_back(_tmp);
      }

      if (line_by_space.size() > 1) {
        for (size_t i = 0; i < line_by_space.size() - 1; ++i)
          adapter_name += line_by_space[i] + " ";
        adapter_seq = line_by_space.back();

        if (adapter_seq.size() > 32)
          throw runtime_error("adapter too long. Maximum adapter size is 32bp: "
                              + adapter_seq);
        adapter_hash = 0;
        char c;
        for (size_t i = 0; i < adapter_seq.size(); ++i) {
          c = adapter_seq[i];
          if (c != 'A' && c != 'C' && c != 'T' && c != 'G')
            throw runtime_error("Bad adapter (non-ATGC characters): "
                                + adapter_seq);

          adapter_hash = (adapter_hash << 2) | actg_to_2bit(c);
        }
      }

      // store information
      adapter_names.push_back(adapter_name);
      adapter_seqs.push_back(adapter_seq);
      adapter_hashes.push_back(adapter_hash);
      line_by_space.clear();
    }
  }
  in.close();
}

void
FalcoConfig::read_contaminants_file() {

  ifstream in(contaminants_file);
  if (!in)
    throw runtime_error("contaminants file not found: " + contaminants_file);

  vector<string> line_by_space;

  // The contaminants file has a space separated name, and the last
  // instance is the biological sequence
  string line;
  while (getline(in, line)) {
    if (is_content_line(line)) {
      string contaminant_seq = "";
      istringstream iss(line);
      string token;
      while (iss >> token)
        line_by_space.push_back(token);

      if (line_by_space.size() > 1) {
        string contaminant_name;
        for (size_t i = 0; i < line_by_space.size() - 1; ++i)
          contaminant_name += line_by_space[i] + " ";
        const string contaminent_seq(line_by_space.back());
        contaminants.push_back(make_pair(contaminant_name, contaminant_seq));
      }
      line_by_space.clear();
    }
  }
}

