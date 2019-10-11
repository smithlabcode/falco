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

#include "smithlab_utils.hpp"
#include "OptionParser.hpp"
#include "Config.hpp"
#include "FastqStats.hpp"
#include "StreamReader.hpp"
#include "HtmlMaker.hpp"
#include <fstream>
#include <chrono>

using std::string;
using std::runtime_error;
using std::cerr;
using std::endl;
using std::vector;
using std::ofstream;
using std::to_string;

using std::chrono::system_clock;

// Function to print progress nicely
static inline void
log_process(const std::string &s) {
  auto tmp = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::string time_fmt = std::string(std::ctime(&tmp));
  time_fmt.pop_back();
  std::cerr << "[" << time_fmt << "] " << s << "\n";
}

// Read any file type until the end and logs progress
template <typename T> void
read_stream_into_stats(T &in,
                       FastqStats &stats,
                       Config &config) {
  in.load();

  // Fast way to report # of reads without modulo arithmetics
  const size_t num_reads_to_log = 1000000;
  size_t next_read = num_reads_to_log;

  // Read record by record
  while (in >> stats) {
    if (!config.quiet) {
      if (stats.num_reads == next_read) {
        log_process ("Processed " +
                    to_string(stats.num_reads / num_reads_to_log) +
                    "M reads");
        next_read += num_reads_to_log;
      }
    }
  }

  // if I could not get tile information from read names, I need to tell this to
  // config so it does not output tile data on the summary or html
  if (in.tile_ignore) {
    config.do_tile = false;
  }
}

int main(int argc, const char **argv) {
  /****************** COMMAND LINE OPTIONS ********************/
  bool help = false;
  bool version = false;

  // skip outputs
  bool skip_text = false;
  bool skip_html = false;

  Config config;
  string forced_file_format;
  string outdir = "";

  OptionParser opt_parse(argv[0],
                         "A high throughput sequence QC analysis tool",
                         "seqfile1 seqfile2 ... seqfileN");

  opt_parse.add_opt("-help", 'h', "print this help file adn exit",
                     false, help);

  opt_parse.add_opt("-version", 'v', "print the program version and exit",
                     false, version);

  opt_parse.add_opt("-outdir", 'o', "Create all output files in the specified "
                    "output directory. If not provided, files will be created "
                    "in the same directory as the input file.",
                    false, outdir);

  opt_parse.add_opt("-casava", 'C',
                    "Files come from raw casava output (currently ignored)",
                    false, config.casava);

  opt_parse.add_opt("-nano", 'n',
                    "Files come from fast5 nanopore sequences",
                    false, config.nanopore);

  opt_parse.add_opt("-nofilter", 'F',
                    "If running with --casava do not remove poor quality "
                    "sequences (currently ignored)",
                    false, config.nofilter);

  opt_parse.add_opt("-noextract", 'e',
                    "If running with --casava do not remove poor quality "
                    " sequences (currently ignored)",
                    false, config.casava);

  opt_parse.add_opt("-nogroup", 'g',
                    "Disable grouping of bases for reads >50bp",
                    false, config.nogroup);

  opt_parse.add_opt("-format", 'f',
                    "Force file format",
                    false, forced_file_format);

  opt_parse.add_opt("-threads", 't',
                    "Specifies number of simultaneous files ",
                    false, config.threads);

  opt_parse.add_opt("-contaminants", 'c',
                    "Non-default filer with a list of contaminants",
                    false, config.contaminants_file);

  opt_parse.add_opt("-adapters", 'a',
                    "Non-default file with a list of adapters",
                    false, config.contaminants_file);

  opt_parse.add_opt("-limits", 'l',
                    "Non-default file with limits and warn/fail criteria",
                    false, config.contaminants_file);

  opt_parse.add_opt("-kmer", 'k',
                    "k-mer size (default = 7, max = 10)", false,
                    config.kmer_size);

  opt_parse.add_opt("-skip-text", 'T',
                    "Skip generating text file (Default = false)", false,
                    skip_text);

  opt_parse.add_opt("-skip-html", 'H',
                    "Skip generating HTML file (Default = false)", false,
                    skip_html);

  opt_parse.add_opt("-quiet", 'q', "print more run info", false, config.quiet);
  opt_parse.add_opt("-dir", 'd', "directory in which to create temp files", 
                     false, config.quiet);

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

  for (auto file : leftover_args) {
    auto file_start = system_clock::now();
    config.filename = file;

    // if format was not provided, we have to guess it by the filename
    if (!forced_file_format.empty())
      config.format = forced_file_format;
    else
      config.format = "";

    /****************** BEGIN PROCESSING CONFIG ******************/
    // define file type, read limits, adapters, contaminants and fail early if
    // any of the required files is not present
    config.setup(); 

    /****************** END PROCESSING CONFIG *******************/
    if (!config.quiet)
      log_process("Started reading file " +config.filename);

    // Allocates vectors to summarize data
    FastqStats stats = FastqStats();

    // Initializes a reader given the file format
    if (config.is_sam) {
      log_process("reading file as sam format");
      SamReader in(config, stats.kNumBases);
      read_stream_into_stats(in,stats,config);
    }
#ifdef USE_HTS
    else if (config.is_bam) {
      log_process("reading file as bam format");
      BamReader in (config, stats.kNumBases);
      read_stream_into_stats(in,stats,config);
    }
#endif

#ifdef USE_ZLIB
    else if (config.is_fastq_gz) {
      log_process("reading file as gzipped fastq format");
      GzFastqReader in (config, stats.kNumBases);
      read_stream_into_stats(in,stats,config);
    }
#endif
    else if (config.is_fastq) {
      log_process("reading file as uncompressed fastq format");
      FastqReader in(config, stats.kNumBases);
      read_stream_into_stats(in, stats, config);
    }
    else {
      throw runtime_error("Cannot recognize file format for file "
                          + config.filename);
    }

    if (!config.quiet)
      log_process("Finished reading file");

    if (!config.quiet)
      log_process("Summarizing data");

    // This function has to be called before writing to output. This is where we
    // calculate all the summary statistics that will be written to output.
    stats.summarize(config);

    /************************ WRITE TO OUTPUT *****************************/
    if (!skip_text) {
      string outfile = file;
      if (!outdir.empty())
        outfile = outdir + "/" + config.filename_stripped;
      outfile += "_qc_data.txt";

      if (!config.quiet)
        log_process("Writing text data to " + outfile);

      // define output
      ofstream of(outfile);

      // Write
      stats.write(of, config);
    }

    /************************ WRITE TO HTML *****************************/
    if (!skip_html) {
      // Take the html template and populate it with stats data.
      // Use config to know which parts not to write
      HtmlMaker html_maker(config.html_file);
      html_maker.write(stats, config);

      // Decide html filename based on input
      string htmlfile = file;
      if (!outdir.empty())
        htmlfile = outdir + "/" + config.filename_stripped;
      htmlfile += "_report.html";

      if (!config.quiet)
        log_process("Writing HTML report to " + htmlfile);

      // Start writing
      ofstream html(htmlfile);
      html << html_maker.sourcecode;
    }

    /************************** TIME SUMMARY *****************************/
    if (!config.quiet) {
      cerr << "Elapsed time for file " << file << ": " <<
      std::chrono::duration_cast<std::chrono::seconds>(
          system_clock::now() - file_start).count()
      << "s\n";
    }
  }

  return EXIT_SUCCESS;
}

