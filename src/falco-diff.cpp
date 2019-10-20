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

#include <fstream>
#include <chrono>

#include "smithlab_utils.hpp"
#include "OptionParser.hpp"
#include "FalcoConfig.hpp"
#include "FastqStats.hpp"
#include "StreamReader.hpp"
#include "HtmlMaker.hpp"

using std::string;
using std::runtime_error;
using std::cerr;
using std::endl;
using std::vector;
using std::ofstream;
using std::to_string;

using std::chrono::system_clock;
using std::chrono::duration_cast;
using time_point = std::chrono::time_point<std::chrono::system_clock>;

// Function to get seconds elapsed in program
static size_t
get_seconds_since(const time_point &start_time) {
  auto current_time = system_clock::now();
  auto time_difference = current_time - start_time;
  return duration_cast<std::chrono::seconds>(time_difference).count();
}

// Function to print progress nicely
static inline void
log_process(const string &s) {
  auto tmp = system_clock::to_time_t(system_clock::now());
  string time_fmt(std::ctime(&tmp));
  time_fmt.pop_back();
  cerr << "[" << time_fmt << "] " << s << endl;
}

// Function to check existance of directory
static bool dir_exists(const string &path) {
  struct stat info;
  if (stat(path.c_str(), &info ) != 0)
    return false;
  else if (info.st_mode & S_IFDIR)
    return true;
  return false;
}

// Read any file type until the end and logs progress
template <typename T> void
read_stream_into_stats(T &in, FastqStats &stats, FalcoConfig &falco_config) {

  in.load();

  // fast way to report # of reads without modulo arithmetics
  const size_t num_reads_to_log = 1000000; // million
  size_t next_read = num_reads_to_log;

  // Read record by record
  while (in >> stats) {
    if (!falco_config.quiet) {
      if (stats.num_reads == next_read) {
        log_process("Processed " +
                    to_string(stats.num_reads/num_reads_to_log) +
                    "M reads");
        next_read += num_reads_to_log;
      }
    }
  }

  // if I could not get tile information from read names, I need to tell this to
  // config so it does not output tile data on the summary or html
  if (in.tile_ignore)
    falco_config.do_tile = false;
}

int main(int argc, const char **argv) {

  try {
    bool help = false;
    bool version = false;

    // skip outputs
    bool skip_text = false;
    bool skip_html = false;
    bool skip_short_summary = false;

    FalcoConfig falco_config;
    string forced_file_format;

    string outdir;
    const string outdir_description =
      "Create all output files in the specified output directory. If not"
      "provided, files will be created in the same directory as the input "
      "file.";

    static const string description =
      "A high throughput sequence QC analysis tool";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], description, "<seqfile1> <seqfile2> ...");
    opt_parse.add_opt("-help", 'h', "print this help file and exit", false,
                        help);
    opt_parse.add_opt("-version", 'v', "print the program version and exit",
                      false, version);
    opt_parse.add_opt("-outdir", 'o', outdir_description, false, outdir);
    opt_parse.add_opt("-format", 'f', "Force file format", false,
                      forced_file_format);
    opt_parse.add_opt("-threads", 't', "Specifies number of simultaneous files",
                      false, falco_config.threads);
    opt_parse.add_opt("-contaminants", 'c',
                      "Non-default filer with a list of contaminants",
                      false, falco_config.contaminants_file);
    opt_parse.add_opt("-adapters", 'a',
                      "Non-default file with a list of adapters",
                      false, falco_config.contaminants_file);
    opt_parse.add_opt("-limits", 'l',
                      "Non-default file with limits and warn/fail criteria",
                      false, falco_config.contaminants_file);
    opt_parse.add_opt("-kmer", 'k', "k-mer size (default = 7, max = 10)", false,
                      falco_config.kmer_size);
    opt_parse.add_opt("-skip-text", 'T', "Skip generating text file "
                      "(Default = false)", false, skip_text);
    opt_parse.add_opt("-skip-html", 'H', "Skip generating HTML file "
                      "(Default = false)", false, skip_html);
    opt_parse.add_opt("-skip-short-summary", 'S', "Skip short summary"
                      "(Default = false)", false, skip_short_summary);
    opt_parse.add_opt("-quiet", 'q', "print more run info", false, falco_config.quiet);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      cerr << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      cerr << opt_parse.option_missing_message() << endl;
      return EXIT_SUCCESS;
    }

    if (leftover_args.size() != 2) {
      cerr << "Number of arguments should be two!\n";
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }

    if (!outdir.empty() && !dir_exists(outdir)) {
      cerr << "output directory does not exist: " << outdir << endl;
      return EXIT_FAILURE;
    }
    const string file1 = leftover_args[0],
                 file2 = leftover_args[1];
    /****************** END COMMAND LINE OPTIONS ********************/
    // Test if both summary files can be opened
    ifstream fa(file1);
    if (!fa)
      throw runtime_error("File not found: " + file1);

    ifstream fb(file2);
    if (!fb)
      throw runtime_error("File not found: " + file2);

    // Make FastqStats objects from files
    FastqStats sa, sb;
    sa.read(file1);
    sb.read(file2);
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


