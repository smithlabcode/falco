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

#include "FalcoConfig.hpp"
#include "FastqStats.hpp"
#include "HtmlMaker.hpp"
#include "Module.hpp"
#include "OptionParser.hpp"
#include "smithlab_utils.hpp"

#include <chrono>
#include <fstream>
#include <sys/stat.h>

// Function to get seconds elapsed in program
static size_t
get_seconds_since(
  const std::chrono::time_point<std::chrono::system_clock> &start_time) {
  auto current_time = std::chrono::system_clock::now();
  auto time_difference = current_time - start_time;
  return std::chrono::duration_cast<std::chrono::seconds>(time_difference)
    .count();
}

// Function to print progress nicely
static inline void
log_process(const std::string &s) {
  auto tmp = system_clock::to_time_t(system_clock::now());
  std::string time_fmt(std::ctime(&tmp));
  time_fmt.pop_back();
  std::cerr << "[" << time_fmt << "] " << s << '\n';
}

// Function to check existance of directory
static bool
dir_exists(const std::string &path) {
  struct stat info;
  if (stat(path.c_str(), &info) != 0)
    return false;
  else if (info.st_mode & S_IFDIR)
    return true;
  return false;
}

// ***********************************************************
// ************ FUNCTIONS TO PARSE TEXT OUTPUT ***************
// ***********************************************************

// skip lines that are empty or comments
inline bool
is_skippable(const std::string &line) {
  return (line.empty() || (line[0] == '#'));
}

inline bool
is_module_start(const std::string &line) {
  return (line.rfind(">>", 0) == 0);
}

inline bool
is_module_end(const std::string &line) {
  return (line == ">>END_MODULE");
}

std::string
get_module_name(const std::string &line) {
  if (line.size() < 2 || line[0] != '>' || line[1] != '>')
    throw std::runtime_error("Bad module name line: " + line);

  // two things: take the >> and the "pass/warn/fail" at the end
  std::string ans = line.substr(2);
  size_t last_tab = ans.find_last_of("\t");
  ans = ans.substr(0, last_tab);

  return ans;
}

// prints in tsv the value, a/b and
void
print_measure_and_diff(std::ostream &os, const std::string measure, double v1,
                       double v2) {
  os << measure << "\t" << v1 << ":" << v2 << "\t" << abs(v1 - v2) << "\n";
}

template <typename T>
void
parse_module(T m1, T m2, std::istream &lhs, std::istream &rhs) {
  std::string line_lhs, line_rhs;
  while (getline(lhs, line_lhs) && getline(rhs, line_rhs) &&
         !is_module_end(line_lhs) && !is_module_end(line_rhs)) {

    if (!is_skippable(line_lhs) && !is_skippable(line_rhs)) {
      m1.read_data_line(line_lhs);
      m2.read_data_line(line_rhs);
    }
  }
}

void
compare_fastqcs(std::istream &lhs, std::istream &rhs, std::ostream &os,
                const FalcoConfig &config) {
  std::string line_lhs, line_rhs, module_name;
  while (!lhs.eof() && !rhs.eof()) {
    getline(lhs, line_lhs);
    getline(rhs, line_rhs);

    if (!is_skippable(line_lhs)) {
      if (is_skippable(line_rhs))
        throw std::runtime_error("lines not aligned (rhs skippable, lhs not");

      if (is_module_start(line_lhs)) {
        if (line_lhs != line_rhs)
          throw std::runtime_error("modules not aligned. lhs = " + line_lhs +
                                   "rhs = " + line_rhs);

        std::string module_name = get_module_name(line_lhs);
        // Basic statistics
        if (module_name == ModuleBasicStatistics::module_name) {
          // define
          ModuleBasicStatistics m1(config), m2(config);
          // read
          std::cerr << "parsing basic statistics\n";
          parse_module(m1, m2, lhs, rhs);

          // write
          os << "#" << ModuleBasicStatistics::module_name << "\n";
          os << "LHS file:\t" << m1.filename_stripped << "\n";
          os << "RHS file:\t" << m2.filename_stripped << "\n";
          os << "#measure\ta:b\tabs_diff\n";
          print_measure_and_diff(os, "total_sequences", m1.total_sequences,
                                 m2.total_sequences);
          print_measure_and_diff(os, "min_read_length", m1.min_read_length,
                                 m2.min_read_length);
          print_measure_and_diff(os, "max_read_length", m1.max_read_length,
                                 m2.max_read_length);
          print_measure_and_diff(os, "gc_percent", m1.avg_gc, m2.avg_gc);
        }

        // Per base sequence quality
        if (module_name == ModulePerBaseSequenceQuality::module_name) {
          ModulePerBaseSequenceQuality m1(config), m2(config);
          std::cerr << "parsing per base sequence quality\n";
          // read
          while (getline(lhs, line_lhs) && getline(rhs, line_rhs) &&
                 !is_module_end(line_lhs) && !is_module_end(line_rhs)) {

            if (!is_skippable(line_lhs) && !is_skippable(line_rhs)) {
              m1.read_data_line(line_lhs);
              m2.read_data_line(line_rhs);
            }
          }
        }
      }
    }

    // check if line in other file is also a comment/skippable
    else if (!is_skippable(line_rhs)) {
      throw std::runtime_error("lines not aligned (lhs skippable, rhs not");
    }
  }

  if (!lhs.eof() || !rhs.eof())
    throw std::runtime_error(
      "the two files do not have the same number of lines!");
}

int
main(int argc, const char **argv) {
  try {
    const time_point file_start_time = system_clock::now();
    bool help = false;
    bool VERBOSE = false;

    // skip outputs
    bool skip_text = false;
    bool skip_html = false;
    bool skip_short_summary = false;

    FalcoConfig config;

    std::string outdir;
    const std::string outdir_description =
      "Create all output files in the specified output directory. If not"
      "provided, files will be created in the same directory as the input "
      "file.";

    static const std::string description =
      "A tool for quality comparison of two samples";

    /****************** COMMAND LINE OPTIONS ********************/
    OptionParser opt_parse(argv[0], description,
                           "<fastqc_data_1.txt> <fastqc_data_2.txt>");
    opt_parse.add_opt("-help", 'h', "print this help file and exit", false,
                      help);
    opt_parse.add_opt("-outdir", 'o', outdir_description, false, outdir);
    opt_parse.add_opt("-skip-text", 'T',
                      "Skip generating text file "
                      "(Default = false)",
                      false, skip_text);
    opt_parse.add_opt("-skip-html", 'H',
                      "Skip generating HTML file "
                      "(Default = false)",
                      false, skip_html);
    opt_parse.add_opt("-skip-short-summary", 'S',
                      "Skip short summary"
                      "(Default = false)",
                      false, skip_short_summary);
    opt_parse.add_opt("-verbose", 'v', "print more run info", false, VERBOSE);

    std::vector<std::string> leftover_args;
    opt_parse.parse(argc, argv, leftover_args);
    if (argc == 1 || opt_parse.help_requested()) {
      std::cerr << opt_parse.help_message() << '\n'
                << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.about_requested()) {
      std::cerr << opt_parse.about_message() << '\n';
      return EXIT_SUCCESS;
    }
    if (opt_parse.option_missing()) {
      std::cerr << opt_parse.option_missing_message() << '\n';
      return EXIT_SUCCESS;
    }

    if (leftover_args.size() != 2) {
      std::cerr << "Number of arguments should be two!\n";
      std::cerr << opt_parse.help_message() << '\n';
      return EXIT_SUCCESS;
    }

    if (!outdir.empty() && !dir_exists(outdir)) {
      std::cerr << "output directory does not exist: " << outdir << '\n';
      return EXIT_FAILURE;
    }
    const std::string file1 = leftover_args[0];
    const std::string file2 = leftover_args[1];

    // read the two files and check if they exist
    std::ifstream lhs(file1), rhs(file2);
    if (!lhs)
      throw std::runtime_error("File not found: " + file1);
    if (!rhs)
      throw std::runtime_error("File not found: " + file2);

    /****************** END COMMAND LINE OPTIONS ********************/

    compare_fastqcs(lhs, rhs, std::cout, config);
    if (VERBOSE)
      std::cerr << "Elapsed time: " << get_seconds_since(file_start_time)
                << "s\n";

    return EXIT_SUCCESS;
  }
  catch (const std::runtime_error &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
