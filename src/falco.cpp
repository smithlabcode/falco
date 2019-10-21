/* falco: quality control for sequencing read files
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
#include "Module.hpp"

using std::string;
using std::runtime_error;
using std::cerr;
using std::endl;
using std::vector;
using std::ofstream;
using std::ostream;
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

// Function to check existance of directort size_t last_slash_idx =
// filename.rfind('\\');
static bool dir_exists(const string &path) {
  struct stat info;
  if (stat(path.c_str(), &info ) != 0)
    return false;
  else if (info.st_mode & S_IFDIR)
    return true;
  return false;
}

// Read any file type until the end and logs progress
// in is an StreamReader object type
template <typename T> void
read_stream_into_stats(T &in, FastqStats &stats, FalcoConfig &falco_config) {
  // open file
  in.load();

  // fast way to report # of reads without modulo arithmetics
  const size_t num_reads_to_log = 1000000; // million
  size_t next_read = num_reads_to_log;

  // Read record by record
  const bool quiet = falco_config.quiet;
  while (in >> stats) {
    if (!quiet) {
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
  if (in.tile_ignore) {
    falco_config.do_tile = false;
  }
}

// Write module content into html maker if requested
template <typename T> void
write_if_requested(T module,
                   const FastqStats &stats,
                   bool requested,
                   bool skip_text,
                   bool skip_html,
                   bool skip_short_summary,
                   const string &filename,
                   ostream &summary_txt,
                   ostream &qc_data_txt,
                   HtmlMaker &html_maker) {

  html_maker.put_comment(module.placeholder_cs, module.placeholder_ce,
                         requested);

  // If module has not been requested we put nothing where the data is
  if (!requested) {
    // puts the actual data (table, graph, etc)
    html_maker.put_data(module.placeholder_data,
                        "");
    return;
  }


  // calculates module summary, mandatory before writing
  module.summarize(stats);

  // writes what was requested
  if (!skip_short_summary) module.write_short_summary(summary_txt, filename);
  if (!skip_text) module.write(qc_data_txt);
  if (!skip_html) {
    // puts the module name 
    html_maker.put_data(module.placeholder_name, module.module_name);

    // puts the grade
    html_maker.put_data(module.placeholder_grade, module.grade);

    // puts the actual data (table, graph, etc)
    html_maker.put_data(module.placeholder_data,
                        module.html_data);
  }
}

void
write_results(const FalcoConfig &falco_config,
              const FastqStats &stats,
              bool skip_text,
              bool skip_html,
              bool skip_short_summary,
              const string &outdir) {

  // Here we open the short summary ofstream
  ofstream summary_txt;
  if (!skip_short_summary) {
    string summary_file = outdir + "/summary.txt";
    summary_txt.open(summary_file.c_str(), std::ofstream::binary);
  }

  // Here we open the full text summary
  ofstream qc_data_txt;
  if (!skip_text) {
    string qc_data_file = falco_config.filename;
    qc_data_file = outdir + "/fastqc_data.txt";
    qc_data_txt.open(qc_data_file.c_str(), std::ofstream::binary);
    // put header
    qc_data_txt << "##Falco\t0.11.8\n";
  }

  // Here we open the html ostream and maker object
  HtmlMaker html_maker(falco_config.html_file);
  ofstream html;
  if (!skip_html) {
    // Decide html filename based on input
    string html_file = outdir + "/fastqc_report.html";

    if (!falco_config.quiet)
      log_process("Writing HTML report to " + html_file);

    html_maker.put_file_details(falco_config);
    html.open(html_file.c_str(), std::ofstream::binary);
  }

  // Now we create modules if requested, summarize them and kill them
  //  Basic Statistics
  write_if_requested(ModuleBasicStatistics(falco_config),
                     stats,
                     true,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);


  //  Per base sequence quality
  write_if_requested(ModulePerBaseSequenceQuality(falco_config),
                     stats,
                     falco_config.do_quality_sequence,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);

  //  Per tile sequence quality
  write_if_requested(ModulePerTileSequenceQuality(falco_config),
                     stats,
                     falco_config.do_tile,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);

  //  Per sequence quality scores
  write_if_requested(ModulePerSequenceQualityScores(falco_config),
                     stats,
                     falco_config.do_quality_sequence,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);

  //  Per base sequence content
  write_if_requested(ModulePerBaseSequenceContent(falco_config),
                     stats,
                     falco_config.do_sequence,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);
  //  Per sequence GC content
  write_if_requested(ModulePerSequenceGCContent(falco_config),
                     stats,
                     falco_config.do_gc_sequence,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);

  //  Per base N content
  write_if_requested(ModulePerBaseNContent(falco_config),
                     stats,
                     falco_config.do_n_content,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);

  //  Sequence Length Distribution
  write_if_requested(ModuleSequenceLengthDistribution(falco_config),
                     stats,
                     falco_config.do_sequence_length,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);

  //  Sequence Duplication Levels
  write_if_requested(ModuleSequenceDuplicationLevels(falco_config),
                     stats,
                     falco_config.do_duplication,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);

  //  Overrepresented sequences
  write_if_requested(ModuleOverrepresentedSequences(falco_config),
                     stats,
                     falco_config.do_overrepresented,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);
  //  Adapter Content
  write_if_requested(ModuleAdapterContent(falco_config),
                     stats,
                     falco_config.do_adapter,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);
  //  Kmer Content
  write_if_requested(ModuleKmerContent(falco_config),
                     stats,
                     falco_config.do_kmer,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped,
                     summary_txt, qc_data_txt,
                     html_maker);


  if (!skip_html)
    html << html_maker.html_boilerplate;
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

    string outdir, tmpdir;
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
    opt_parse.add_opt("-casava", 'C',
                      "Files come from raw casava output (currently ignored)",
                      false, falco_config.casava);
    opt_parse.add_opt("-nano", 'n', "Files come from fast5 nanopore sequences",
                      false, falco_config.nanopore);
    opt_parse.add_opt("-nofilter", 'F', "If running with --casava do not "
                      "sequences (currently ignored)", false,
                      falco_config.nofilter);
    opt_parse.add_opt("-noextract", 'e', "If running with --casava do not "
                      "remove poor quality sequences (currently ignored)",
                      false, falco_config.casava);
    opt_parse.add_opt("-nogroup", 'g', "Disable grouping of bases for "
                      "reads >50bp", false, falco_config.nogroup);
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
    opt_parse.add_opt("-dir", 'd', "directory in which to create temp files",
                      false, tmpdir);

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

    if (leftover_args.size() == 0) {
      cerr << opt_parse.help_message() << endl;
      return EXIT_SUCCESS;
    }

    if (!outdir.empty() && !dir_exists(outdir)) {
      cerr << "output directory does not exist: " << outdir << endl;
      return EXIT_FAILURE;
    }
    const vector<string> all_seq_filenames(leftover_args);

    /****************** END COMMAND LINE OPTIONS ********************/
    for (auto filename : all_seq_filenames) {

      const time_point file_start_time = system_clock::now();

      falco_config.filename = filename;

      // if format was not provided, we have to guess it by the filename
      if (!forced_file_format.empty())
        falco_config.format = forced_file_format;
      else
        falco_config.format = "";

      /****************** BEGIN PROCESSING CONFIG ******************/
      // define file type, read limits, adapters, contaminants and fail
      // early if any of the required files is not present
      //
      falco_config.setup();

      /****************** END PROCESSING CONFIG *******************/
      if (!falco_config.quiet)
        log_process("Started reading file " + falco_config.filename);
      FastqStats stats; // allocate all space to summarize data

      // Initializes a reader given the file format
      if (falco_config.is_sam) {
        log_process("reading file as sam format");
        SamReader in(falco_config, stats.kNumBases);
        read_stream_into_stats(in, stats, falco_config);
      }
#ifdef USE_HTS
      else if (falco_config.is_bam) {
        log_process("reading file as bam format");
        BamReader in(falco_config, stats.kNumBases);
        read_stream_into_stats(in, stats, falco_config);
      }
#endif

      else if (falco_config.is_fastq_gz) {
        log_process("reading file as gzipped fastq format");
        GzFastqReader in(falco_config, stats.kNumBases);
        read_stream_into_stats(in,stats,falco_config);
      }
      else if (falco_config.is_fastq) {
        log_process("reading file as uncompressed fastq format");
        FastqReader in(falco_config, stats.kNumBases);
        read_stream_into_stats(in, stats, falco_config);
      }
      else {
        throw runtime_error("Cannot recognize file format for file "
                            + falco_config.filename);
      }

      if (!falco_config.quiet) {
        log_process("Finished reading file");
      }

      stats.summarize();

      // if oudir is empty we will set it as the filename path
      string cur_outdir;
      if (outdir.empty()) {
        const size_t last_slash_idx = filename.rfind('/');
        // if file was given with relative path in the current dir, we set a dot
        if (last_slash_idx == string::npos) {
          cur_outdir = ".";
        } else {
          cur_outdir = falco_config.filename.substr(0, last_slash_idx);
        }
      } else {
        cur_outdir = outdir;
      }

      // Write results
      log_process("Writing results");
      write_results(falco_config, stats, skip_text, skip_html,
                   skip_short_summary, cur_outdir);

      /************************** TIME SUMMARY *****************************/
      if (!falco_config.quiet)
        cerr << "Elapsed time for file " << filename << ": "
             << get_seconds_since(file_start_time) << "s" << endl;
    }
  }
  catch (const runtime_error &e) {
    cerr << e.what() << endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}


