/* falco: quality control for sequencing read files
 *
 * Copyright (C) 2019-2022 Guilherme De Sena Brandine and
 *                         Andrew D. Smith
 * Authors: Guilherme De Sena Brandine, Andrew Smith, and Masaru Nakajima
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

#include <chrono>
#include <filesystem>
#include <fstream>

#include "FalcoConfig.hpp"
#include "FastqStats.hpp"
#include "HtmlMaker.hpp"
#include "Module.hpp"
#include "OptionParser.hpp"
#include "StreamReader.hpp"
#include "smithlab_utils.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::runtime_error;
using std::string;
using std::to_string;
using std::vector;

namespace fs = std::filesystem;

using std::chrono::duration_cast;
using std::chrono::system_clock;
using time_point = std::chrono::time_point<std::chrono::system_clock>;

// Function to get seconds elapsed in program
static size_t
get_seconds_since(const time_point &start_time) {
  const auto current_time = system_clock::now();
  const auto time_difference = current_time - start_time;
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
static bool
dir_exists(const string &path) {
  return fs::exists(path) && fs::is_directory(path);
}

// Read any file type until the end and logs progress
// in is an StreamReader object type
template <typename T>
void
read_stream_into_stats(T &in, FastqStats &stats, FalcoConfig &falco_config) {
  // open file
  size_t file_size = std::max(in.load(), static_cast<size_t>(1));
  size_t tot_bytes_read = 0;

  // Read record by record
  const bool quiet = falco_config.quiet;
  ProgressBar progress(file_size, "running falco");
  if (!quiet)
    progress.report(cerr, 0);
  while (in.read_entry(stats, tot_bytes_read)) {
    if (!quiet && progress.time_to_report(tot_bytes_read))
      progress.report(cerr, tot_bytes_read);
  }

  // if I could not get tile information from read names, I need to tell this to
  // config so it does not output tile data on the summary or html
  if (in.tile_ignore)
    falco_config.do_tile = false;

  if (!quiet && tot_bytes_read < file_size)
    progress.report(cerr, file_size);
}

// Write module content into html maker if requested
template <typename T>
void
write_if_requested(T module, FastqStats &stats, const bool requested,
                   const bool skip_text, const bool skip_html,
                   const bool skip_short_summary, const string &filename,
                   ostream &summary_txt, ostream &qc_data_txt,
                   HtmlMaker &html_maker) {
  html_maker.put_comment(module.placeholder_cs, module.placeholder_ce,
                         requested);

  // If module has not been requested we put nothing where the data is
  if (!requested) {
    // puts the actual data (table, graph, etc)
    html_maker.put_data(module.placeholder_data, "");
    return;
  }

  // calculates module summary, mandatory before writing
  module.summarize(stats);

  // writes what was requested
  if (!skip_short_summary)
    module.write_short_summary(summary_txt, filename);
  if (!skip_text)
    module.write(qc_data_txt);
  if (!skip_html) {
    // puts the module name
    html_maker.put_data(module.placeholder_name, T::module_name);

    // puts the grade
    html_maker.put_data(module.placeholder_grade, module.grade);

    // puts the actual data (table, graph, etc)
    html_maker.put_data(module.placeholder_data, module.html_data);
  }
}

void
write_results(const FalcoConfig &falco_config, FastqStats &stats,
              const bool skip_text, const bool skip_html,
              const bool skip_short_summary, const bool do_call,
              const string &file_prefix, const string &outdir,
              const string &summary_filename, const string &data_filename,
              const string &report_filename) {

  // Here we open the short summary ofstream
  ofstream summary_txt;
  if (!skip_short_summary) {
    const string summary_file =
      (summary_filename.empty() ? (outdir + "/" + file_prefix + "summary.txt")
                                : (summary_filename));
    summary_txt.open(summary_file.c_str(), std::ofstream::binary);

    if (!summary_txt.good())
      throw runtime_error("Failed to create output summary file: " +
                          summary_file);

    if (!falco_config.quiet)
      log_process("Writing summary to " + summary_file);
  }

  // Here we open the full text summary
  ofstream qc_data_txt;
  if (!skip_text) {
    const string qc_data_file =
      (data_filename.empty() ? (outdir + "/" + file_prefix + "fastqc_data.txt")
                             : (data_filename));
    qc_data_txt.open(qc_data_file.c_str(), std::ofstream::binary);

    if (!qc_data_txt.good())
      throw runtime_error("Failed to create output data file: " + qc_data_file);

    if (!falco_config.quiet)
      log_process("Writing text report to " + qc_data_file);

    // put header
    qc_data_txt << "##Falco\t" + FalcoConfig::FalcoVersion + "\n";
    if (do_call)
      qc_data_txt << "##Call\t" << falco_config.call << "\n";
  }

  // Here we open the html ostream and maker object
  HtmlMaker html_maker = HtmlMaker();
  ofstream html;
  if (!skip_html) {
    // Decide html filename based on input
    const string html_file =
      (report_filename.empty()
         ? (outdir + "/" + file_prefix + "fastqc_report.html")
         : (report_filename));

    if (!falco_config.quiet)
      log_process("Writing HTML report to " + html_file);

    html.open(html_file.c_str(), std::ofstream::binary);
    if (!html.good())
      throw runtime_error("Failed to create output HTML report file: " +
                          html_file);

    html_maker.put_file_details(falco_config);
  }

  // Now we create modules if requested, summarize them and kill them
  //  Basic Statistics
  write_if_requested(ModuleBasicStatistics(falco_config), stats, true,
                     skip_text, skip_html, skip_short_summary,
                     falco_config.filename_stripped, summary_txt, qc_data_txt,
                     html_maker);

  //  Per base sequence quality
  write_if_requested(ModulePerBaseSequenceQuality(falco_config), stats,
                     falco_config.do_quality_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per tile sequence quality
  write_if_requested(ModulePerTileSequenceQuality(falco_config), stats,
                     falco_config.do_tile, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per sequence quality scores
  write_if_requested(ModulePerSequenceQualityScores(falco_config), stats,
                     falco_config.do_quality_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per base sequence content
  write_if_requested(ModulePerBaseSequenceContent(falco_config), stats,
                     falco_config.do_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);
  //  Per sequence GC content
  write_if_requested(ModulePerSequenceGCContent(falco_config), stats,
                     falco_config.do_gc_sequence, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Per base N content
  write_if_requested(ModulePerBaseNContent(falco_config), stats,
                     falco_config.do_n_content, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Sequence Length Distribution
  write_if_requested(ModuleSequenceLengthDistribution(falco_config), stats,
                     falco_config.do_sequence_length, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Sequence Duplication Levels
  write_if_requested(ModuleSequenceDuplicationLevels(falco_config), stats,
                     falco_config.do_duplication, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  //  Overrepresented sequences
  write_if_requested(ModuleOverrepresentedSequences(falco_config), stats,
                     falco_config.do_overrepresented, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);
  //  Adapter Content
  write_if_requested(ModuleAdapterContent(falco_config), stats,
                     falco_config.do_adapter, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);
  //  Kmer Content
  write_if_requested(ModuleKmerContent(falco_config), stats,
                     falco_config.do_kmer, skip_text, skip_html,
                     skip_short_summary, falco_config.filename_stripped,
                     summary_txt, qc_data_txt, html_maker);

  if (!skip_html)
    html << html_maker.html_boilerplate;
}

inline bool
file_exists(const string &file_name) {
  return (access(file_name.c_str(), F_OK) == 0);
}

int
main(int argc, const char **argv) {

  try {
    static const string FALCO_VERSION = "falco " + FalcoConfig::FalcoVersion;
    bool help = false;
    bool version = false;

    // skip outputs
    bool skip_text = false;
    bool skip_html = false;
    bool skip_short_summary = false;
    bool do_call = false;
    bool allow_empty_input = false;

    // a tmp boolean to keep compatibility with FastQC
    bool tmp_compatibility_only = false;
    string tmp_compatibility_only_str;

    FalcoConfig falco_config(argc, argv);

    // if defined, read file as the file format specified by the user
    string forced_file_format;

    // output directory in which to put files
    string outdir;

    // custom output filenames (if provided)
    string data_filename = "";
    string report_filename = "";
    string summary_filename = "";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],
                           "A high throughput sequence QC analysis tool",
                           "<seqfile1> <seqfile2> ...");
    opt_parse.set_show_defaults();

    /***************** FASTQC PARAMETERS ********************/
    opt_parse.add_opt("-help", 'h', "Print this help file and exit", false,
                      help);

    opt_parse.add_opt("-version", 'v',
                      "Print the version of the program and exit", false,
                      version);

    opt_parse.add_opt(
      "-outdir", 'o',
      "Create all output files in the specified output directory. "
      "FALCO-SPECIFIC: If the directory does not exists, the "
      "program will create it. "
      "If this option is not set then the output file for each "
      "sequence file is created in the same directory as the "
      "sequence file which was processed.",
      false, outdir);

    opt_parse.add_opt(
      "-casava", '\0',
      "[IGNORED BY FALCO] "
      "Files come from raw casava output. Files in the same sample "
      "group (differing only by the group number) will be analysed "
      "as a set rather than individually. Sequences with the filter "
      "flag set in the header will be excluded from the analysis. "
      "Files must have the same names given to them by casava "
      "(including being gzipped and ending with .gz) otherwise they "
      "won't be grouped together correctly.",
      false, falco_config.casava);

    opt_parse.add_opt(
      "-nano", '\0',
      "[IGNORED BY FALCO] "
      "Files come from nanopore sequences and are in fast5 format. In "
      "this mode you can pass in directories to process and the program "
      "will take in all fast5 files within those directories and produce "
      " a single output file from the sequences found in all files.",
      false, falco_config.nanopore);

    opt_parse.add_opt(
      "-nofilter", '\0',
      "[IGNORED BY FALCO] "
      "If running with --casava then don't remove read flagged by "
      "casava as poor quality when performing the QC analysis. ",
      false, falco_config.nofilter);

    opt_parse.add_opt(
      "-extract", '\0',
      "[ALWAYS ON IN FALCO] "
      "If set then the zipped output file will be uncompressed in "
      "the same directory after it has been created.  By default "
      "this option will be set if fastqc is run in non-interactive "
      "mode.",
      false, tmp_compatibility_only);

    opt_parse.add_opt(
      "-java", 'j',
      "[IGNORED BY FALCO] "
      "Provides the full path to the java binary you want to use to "
      "launch fastqc. If not supplied then java is assumed to be in "
      "your path.",
      false, tmp_compatibility_only_str);

    opt_parse.add_opt(
      "-noextract", '\0',
      "[IGNORED BY FALCO] "
      "Do not uncompress the output file after creating it.  You "
      "should set this option if you do not wish to uncompress "
      "the output when running in non-interactive mode.",
      false, falco_config.casava);

    opt_parse.add_opt(
      "-nogroup", '\0',
      "Disable grouping of bases for reads >50bp. All reports will "
      "show data for every base in the read.  WARNING: When using this "
      "option, your plots may end up a ridiculous size. "
      "You have been warned!",
      false, falco_config.nogroup);

    opt_parse.add_opt(
      "-min_length", '\0',
      "[NOT YET IMPLEMENTED IN FALCO] "
      "Sets an artificial lower limit on the length of the sequence "
      "to be shown in the report.  As long as you set this to a value "
      "greater or equal to your longest read length then this will be "
      "the sequence length used to create your read groups.  This can "
      "be useful for making directly comaparable statistics from "
      "datasets with somewhat variable read lengths. ",
      false, tmp_compatibility_only_str);

    opt_parse.add_opt("-format", 'f',
                      "Bypasses the normal sequence file format detection and "
                      "forces the program to use the specified format.  Valid"
                      "formats are bam, sam, bam_mapped, sam_mapped, fastq, "
                      "fq, fastq.gz or fq.gz.",
                      false, forced_file_format);

    opt_parse.add_opt("-threads", 't',
                      "[NOT YET IMPLEMENTED IN FALCO] "
                      "Specifies the number of files which can be processed "
                      "simultaneously.  Each thread will be allocated 250MB of "
                      "memory so you shouldn't run more threads than your "
                      "available memory will cope with, and not more than "
                      "6 threads on a 32 bit machine",
                      false, falco_config.threads);

    opt_parse.add_opt(
      "-contaminants", 'c',
      "Specifies a non-default file which contains the list of "
      "contaminants to screen overrepresented sequences against. "
      "The file must contain sets of named contaminants in the "
      "form name[tab]sequence.  Lines prefixed with a hash will "
      "be ignored. Default: ",
      false, falco_config.contaminants_file);

    opt_parse.add_opt(
      "-adapters", 'a',
      "Specifies a non-default file which contains the list of "
      "adapter sequences which will be explicity searched against "
      "the library. The file must contain sets of named adapters "
      "in the form name[tab]sequence.  Lines prefixed with a hash "
      "will be ignored. Default: ",
      false, falco_config.adapters_file);

    opt_parse.add_opt(
      "-limits", 'l',
      "Specifies a non-default file which contains a set of criteria "
      "which will be used to determine the warn/error limits for the "
      "various modules.  This file can also be used to selectively "
      "remove some modules from the output all together.  The format "
      "needs to mirror the default limits.txt file found in the "
      "Configuration folder. Default: ",
      false, falco_config.limits_file);

    opt_parse.add_opt(
      "-kmers", 'k',
      "[IGNORED BY FALCO AND ALWAYS SET TO 7] "
      "Specifies the length of Kmer to look for in the Kmer content "
      "module. Specified Kmer length must be between 2 and 10. Default "
      "length is 7 if not specified.",
      false, tmp_compatibility_only_str);

    opt_parse.add_opt(
      "-quiet", 'q',
      "Supress all progress messages on stdout and only report errors.", false,
      falco_config.quiet);

    opt_parse.add_opt(
      "-dir", 'd',
      "[IGNORED: FALCO DOES NOT CREATE TMP FILES] "
      "Selects a directory to be used for temporary files written when "
      "generating report images. Defaults to system temp directory if "
      "not specified. ",
      false, tmp_compatibility_only_str);

    /***************** FALCO ONLY ********************/
    // falco-only options use a single dash as long argument
    // e.g. -skip-text instead of --skip-text. They also
    // use single-characters that do not overlap with the
    // ones used in FastQC, so cannot use:
    // h, v, o, j, f, t, c, a, l, k, q, d

    opt_parse.add_opt("subsample", 's',
                      "[Falco only] makes falco faster "
                      "(but possibly less accurate) by only processing "
                      "reads that are multiple of this value (using 0-based "
                      "indexing to number reads).",
                      false, falco_config.read_step);

    opt_parse.add_opt(
      "bisulfite", 'b',
      "[Falco only] reads are whole genome bisulfite sequencing, and more "
      "Ts and fewer Cs are therefore expected and will be "
      "accounted for in base content.",
      false, falco_config.is_bisulfite);

    opt_parse.add_opt(
      "reverse-complement", 'r',
      "[Falco only] The input is a reverse-complement. All modules will "
      "be tested by swapping A/T and C/G",
      false, falco_config.is_reverse_complement);

    opt_parse.add_opt("skip-data", '\0',
                      "[Falco only] Do not create FastQC data text file.",
                      false, skip_text);

    opt_parse.add_opt("skip-report", '\0',
                      "[Falco only] Do not create FastQC report HTML file.",
                      false, skip_html);

    opt_parse.add_opt("skip-summary", '\0',
                      "[Falco only] Do not create FastQC summary file", false,
                      skip_short_summary);

    opt_parse.add_opt(
      "data-filename", 'D',
      "[Falco only] Specify filename for FastQC data output (TXT). "
      "If not specified, it will be called fastq_data.txt in either "
      "the input file's directory or the one specified in the --output "
      "flag. Only available when running falco with a single input. ",
      false, data_filename);

    opt_parse.add_opt(
      "report-filename", 'R',
      "[Falco only] Specify filename for FastQC report output (HTML). "
      "If not specified, it will be called fastq_report.html in either "
      "the input file's directory or the one specified in the --output "
      "flag. Only available when running falco with a single input.",
      false, report_filename);

    opt_parse.add_opt(
      "summary-filename", 'S',
      "[Falco only] Specify filename for the short summary output (TXT). "
      "If not specified, it will be called fastq_report.html in either "
      "the input file's directory or the one specified in the --output "
      "flag. Only available when running falco with a single input.",
      false, summary_filename);

    opt_parse.add_opt("add-call", 'K',
                      "[Falco only] add the command call call to"
                      " FastQC data output and FastQC report HTML (this"
                      " may break the parse of fastqc_data.txt"
                      " in programs that are very strict about the "
                      " FastQC output format).",
                      false, do_call);

    opt_parse.add_opt(
      "allow-empty-input", '\0',
      "[Falco only] allow empty input files and generate empty output files "
      "without en error state. WARNING: using this option can mask problems in "
      "other parts of a workflow.",
      false, allow_empty_input);

    vector<string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      cerr << opt_parse.help_message() << endl
           << opt_parse.about_message() << endl;
      return EXIT_SUCCESS;
    }
    if (version) {
      cout << FALCO_VERSION << "\n";
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

    if (leftover_args.size() > 1 &&
        (!summary_filename.empty() || !report_filename.empty() ||
         !data_filename.empty())) {
      cerr << "ERROR: summary, report or data filename provided, but"
           << " you are running falco with " << leftover_args.size()
           << " inputs. We cannot allow this because multiple inputs"
           << " require multiple outputs, so they will be resolved by"
           << " the input filenames instead" << endl;
      return EXIT_FAILURE;
    }

    // ADS: make sure all input files are non-empty unless user oks it
    if (!allow_empty_input) {
      for (const auto &fn : leftover_args) {
        std::error_code ec;
        const bool empty_file = std::filesystem::is_empty(fn, ec);
        if (ec) {
          cerr << "Error reading file: " << fn << " (" << ec.message() << ")"
               << endl;
          return EXIT_FAILURE;
        }
        else if (empty_file) {
          cerr << "Input file is empty: " << fn << endl;
          return EXIT_FAILURE;
        }
      }
    }

    if (!outdir.empty()) {
      if (!summary_filename.empty())
        cerr << "[WARNING] specifying custom output directory but also "
             << "custom summary filename. Output " << outdir
             << " will be disregarded and summary file will be saved to "
             << summary_filename << endl;

      if (!data_filename.empty())
        cerr << "[WARNING] specifying custom output directory but also "
             << "custom data filename. Output " << outdir
             << " will be disregarded and data file will be saved to "
             << data_filename << endl;

      if (!report_filename.empty())
        cerr << "[WARNING] specifying custom output directory but also "
             << "custom HTML Report filename. Output " << outdir
             << " will be disregarded and HTML report file will be saved to "
             << report_filename << endl;

      if (!dir_exists(outdir)) {
        if (!falco_config.quiet)
          log_process("creating directory for output: " + outdir);

        // makes directory with r and w permission
        if (mkdir(outdir.c_str(), S_IRUSR | S_IWUSR | S_IXUSR) != 0) {
          cerr << "failed to create directory: " << outdir
               << ". Make sure you have write permissions on it!" << endl;
          return EXIT_FAILURE;
        }
      }
    }
    const vector<string> all_seq_filenames(leftover_args);

    // check if all filenames exist
    bool all_files_exist = true;
    for (size_t i = 0; i < all_seq_filenames.size(); ++i) {
      if (!file_exists(all_seq_filenames[i])) {
        cerr << "ERROR! File does not exist: " << all_seq_filenames[i] << endl;
        all_files_exist = false;
      }
    }

    if (!all_files_exist) {
      throw runtime_error(
        "not all input files exist. Check stderr for detailed list"
        " of files that were not found");
    }

    /****************** END COMMAND LINE OPTIONS ********************/
    for (const auto &filename : all_seq_filenames) {

      const time_point file_start_time = system_clock::now();

      falco_config.filename = filename;

      // if format was not provided, we have to guess it by the filename
      falco_config.format = forced_file_format;

      /****************** BEGIN PROCESSING CONFIG ******************/
      // define file type, read limits, adapters, contaminants and fail
      // early if any of the required files is not present
      //
      falco_config.setup();

      /****************** END PROCESSING CONFIG *******************/
      if (!falco_config.quiet)
        log_process("Started reading file " + falco_config.filename);
      FastqStats stats;  // allocate all space to summarize data

      // Initializes a reader given the file format
      if (falco_config.is_sam) {
        if (!falco_config.quiet)
          log_process("reading file as SAM format");
        SamReader in(falco_config, stats.SHORT_READ_THRESHOLD);
        read_stream_into_stats(in, stats, falco_config);
        stats.adjust_tile_maps_len();
      }
#ifdef USE_HTS
      else if (falco_config.is_bam) {
        if (!falco_config.quiet)
          log_process("reading file as BAM format");
        BamReader in(falco_config, stats.SHORT_READ_THRESHOLD);
        read_stream_into_stats(in, stats, falco_config);
        stats.adjust_tile_maps_len();
      }
#endif

      else if (falco_config.is_fastq_gz) {
        if (!falco_config.quiet)
          log_process("reading file as gzipped FASTQ format");
        GzFastqReader in(falco_config, stats.SHORT_READ_THRESHOLD);
        read_stream_into_stats(in, stats, falco_config);
      }
      else if (falco_config.is_fastq) {
        if (!falco_config.quiet)
          log_process("reading file as uncompressed FASTQ format");
        FastqReader in(falco_config, stats.SHORT_READ_THRESHOLD);
        read_stream_into_stats(in, stats, falco_config);
      }
      else {
        throw runtime_error("Cannot recognize file format for file " +
                            falco_config.filename);
      }

      if (!falco_config.quiet) {
        log_process("Finished reading file");
      }

      stats.summarize();

      // if oudir is empty we will set it as the filename path
      string cur_outdir;

      const size_t last_slash_idx = filename.rfind('/');
      string file_basename = falco_config.filename.substr(last_slash_idx + 1);
      if (outdir.empty()) {
        // if file was given with relative path in the current dir, we set a dot
        if (last_slash_idx == string::npos) {
          cur_outdir = ".";
          file_basename = filename;
        }
        else {
          cur_outdir = falco_config.filename.substr(0, last_slash_idx);
        }
      }
      else {
        cur_outdir = outdir;
      }

      // Write results
      const string file_prefix =
        (all_seq_filenames.size() == 1) ? ("") : (file_basename + "_");
      write_results(falco_config, stats, skip_text, skip_html,
                    skip_short_summary, do_call, file_prefix, cur_outdir,
                    summary_filename, data_filename, report_filename);

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
