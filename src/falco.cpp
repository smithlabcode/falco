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
using std::cout;
using std::endl;
using std::vector;
using std::ifstream;
using std::ofstream;
using std::ostream;
using std::to_string;

using std::chrono::system_clock;
using std::chrono::duration_cast;
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
  size_t file_size = in.load();
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
  if (in.tile_ignore) {
    falco_config.do_tile = false;
  }

  if (tot_bytes_read < file_size && !quiet)
    progress.report(cerr, file_size);

}

// Write module content into html maker if requested
template <typename T> void
write_if_requested(T module,
                   FastqStats &stats,
                   const bool requested,
                   const bool skip_text,
                   const bool skip_html,
                   const bool skip_short_summary,
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
    html_maker.put_data(module.placeholder_name, T::module_name);

    // puts the grade
    html_maker.put_data(module.placeholder_grade, module.grade);

    // puts the actual data (table, graph, etc)
    html_maker.put_data(module.placeholder_data,
                        module.html_data);
  }
}

void
write_results(const FalcoConfig &falco_config,
              FastqStats &stats,
              const bool skip_text,
              const bool skip_html,
              const bool skip_short_summary,
              const bool do_call,
              const string &file_prefix,
              const string &outdir) {

  // Here we open the short summary ofstream
  ofstream summary_txt;
  if (!skip_short_summary) {
    string summary_file = outdir + "/" + file_prefix + "summary.txt";
    summary_txt.open(summary_file.c_str(), std::ofstream::binary);
  }

  // Here we open the full text summary
  ofstream qc_data_txt;
  if (!skip_text) {

    string qc_data_file = falco_config.filename;
    qc_data_file = outdir + "/" + file_prefix + "fastqc_data.txt";
    qc_data_txt.open(qc_data_file.c_str(), std::ofstream::binary);

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
    string html_file = outdir + "/" + file_prefix + "fastqc_report.html";

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

inline bool
file_exists(const string &file_name) {
  return (access(file_name.c_str(), F_OK) == 0);
}

int main(int argc, const char **argv) {

  try {
    static const string FALCO_VERSION = "falco " + FalcoConfig::FalcoVersion;
    bool help = false;
    bool version = false;

    // skip outputs
    bool skip_text = false;
    bool skip_html = false;
    bool skip_short_summary = false;
    bool do_call = false;

    FalcoConfig falco_config(argc, argv);

    // if defined, read file as the file format specified by the user
    string forced_file_format;

    // output directory in which to put files
    string outdir;

    // tmpdir is for fastqc compatibility, not really used
    string tmpdir;

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0],
                              "A high throughput sequence QC analysis tool",
                              "<seqfile1> <seqfile2> ...");
    opt_parse.set_show_defaults();
    opt_parse.add_opt("-help", 'h', "print this help file and exit", false,
                        help);
    opt_parse.add_opt("-version", 'v', "print the program version and exit",
                      false, version);
    opt_parse.add_opt("-outdir", 'o',
      "Create all output files in the specified output directory. If not"
      "provided, files will be created in the same directory as the input "
      "file."
    , false, outdir);

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
    opt_parse.add_opt("-threads", 't', "Specifies number of threads to process "
                                        "simultaneos files in parallel "
                                        "(currently set for compatibility "
                                        "with fastqc. Not yet supported!)",
                      false, falco_config.threads);
    opt_parse.add_opt("-contaminants", 'c',
                      "Non-default filer with a list of contaminants",
                      false, falco_config.contaminants_file);
    opt_parse.add_opt("-adapters", 'a',
                      "Non-default file with a list of adapters",
                      false, falco_config.adapters_file);
    opt_parse.add_opt("-limits", 'l',
                      "Non-default file with limits and warn/fail criteria",
                      false, falco_config.limits_file);
    opt_parse.add_opt("-skip-text", 'T', "Skip generating text file "
                      "(Default = false)", false, skip_text);
    opt_parse.add_opt("-skip-html", 'H', "Skip generating HTML file "
                      "(Default = false)", false, skip_html);
    opt_parse.add_opt("-skip-short-summary", 'S', "Skip short summary"
                      "(Default = false)", false, skip_short_summary);
    opt_parse.add_opt("-quiet", 'q', "print more run info", false,
                      falco_config.quiet);
    opt_parse.add_opt("-dir", 'd', "directory in which to create temp files",
                      false, tmpdir);

    // Falco-specific options
    opt_parse.add_opt("-step", 's', "makes falco faster "
        "(but possibly less accurate) by only processing reads that are multiple "
        "of this value", false, falco_config.read_step);
    opt_parse.add_opt("-bisulfite", 'B',
                      "reads are whole genome bisulfite sequencing, and more "
                      "Ts and fewer Cs are therefore expected and will be "
                      "accounted for in base content (advanced mode)", false,
                      falco_config.is_bisulfite);
    opt_parse.add_opt("-reverse-complement", 'R',
                         "The input is a reverse-complement. All modules will "
                         "be tested by swapping A/T and C/G", false,
                      falco_config.is_reverse_complement);
    opt_parse.add_opt("-add-call", 'K', "add the function call to"
                      " fastqc_data.txt and fastqc-report.html (this"
                      " may break the parse of fastqc_data.txt"
                      " in programs that require rigorous"
                      " FastQC format", false,
                      do_call);
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

    if (!outdir.empty()) {
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
      throw runtime_error("not all input files exist. Check stderr for detailed list"
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
      FastqStats stats; // allocate all space to summarize data

      // Initializes a reader given the file format
      if (falco_config.is_sam) {
        if (!falco_config.quiet)
          log_process("reading file as SAM format");
        SamReader in(falco_config, stats.kNumBases);
        read_stream_into_stats(in, stats, falco_config);
      }
#ifdef USE_HTS
      else if (falco_config.is_bam) {
        if (!falco_config.quiet)
          log_process("reading file as BAM format");
        BamReader in(falco_config, stats.kNumBases);
        read_stream_into_stats(in, stats, falco_config);
      }
#endif

      else if (falco_config.is_fastq_gz) {
        if (!falco_config.quiet)
          log_process("reading file as gzipped FASTQ format");
        GzFastqReader in(falco_config, stats.kNumBases);
        read_stream_into_stats(in,stats,falco_config);
      }
      else if (falco_config.is_fastq) {
        if (!falco_config.quiet)
          log_process("reading file as uncompressed FASTQ format");
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
        }
        else {
          cur_outdir = falco_config.filename.substr(0, last_slash_idx);
        }
      }
      else {
        cur_outdir = outdir;
      }

      // Write results
      const string file_prefix = (all_seq_filenames.size() == 1) ?
                                 ("") : (filename + "_");
      write_results(falco_config, stats, skip_text, skip_html,
                   skip_short_summary, do_call, file_prefix, cur_outdir);

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
