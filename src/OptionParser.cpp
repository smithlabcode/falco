/*
 *    Part of SMITHLAB software
 *
 *    Copyright (C) 2008 Cold Spring Harbor Laboratory,
 *                       University of Southern California and
 *                       Andrew D. Smith
 *
 *    Authors: Andrew D. Smith
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "OptionParser.hpp"

#include <sstream>
#include <cstdlib>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <exception>
#include <cstring>
#include <cctype>
#include <functional>
#include <regex>
#include <iterator>

#include "smithlab_utils.hpp"

using std::vector;
using std::string;
using std::endl;
using std::runtime_error;
using std::regex;
using std::begin;
using std::end;

static const size_t MAX_LINE_LENGTH = 72;

enum {
  SMITHLAB_ARG_INT,    SMITHLAB_ARG_UINT,  SMITHLAB_ARG_LONG,
  SMITHLAB_ARG_ULONG,  SMITHLAB_ARG_FLOAT, SMITHLAB_ARG_DOUBLE,
  SMITHLAB_ARG_STRING, SMITHLAB_ARG_BOOL,  SMITHLAB_ARG_CHAR
};

void
Option::format_option(const string &argument) {
  std::istringstream ss(argument);
  if ((arg_type == SMITHLAB_ARG_INT && !(ss >> *int_value)) ||
      (arg_type == SMITHLAB_ARG_UINT && !(ss >> *uint_value)) ||
      (arg_type == SMITHLAB_ARG_LONG && !(ss >> *long_value)) ||
      (arg_type == SMITHLAB_ARG_ULONG && !(ss >> *ulong_value)) ||
      (arg_type == SMITHLAB_ARG_FLOAT && !(ss >> *float_value)) ||
      (arg_type == SMITHLAB_ARG_DOUBLE && !(ss >> *double_value)) ||
      (arg_type == SMITHLAB_ARG_CHAR && !(ss >> *char_value)))
    throw runtime_error("Invalid argument [" + argument +
                        "] to option [" + format_option_name() + "]");
  else if (arg_type == SMITHLAB_ARG_STRING)
    *string_value = argument;
  else if (arg_type == SMITHLAB_ARG_BOOL) {
    *bool_value = !(*bool_value);
    if (argument == "true" || argument == "on")
      *bool_value = true;
    if (argument == "false" || argument == "off")
      *bool_value = false;
  }
}

using std::numeric_limits;
using std::to_string;

template<class T> string
format_int_like(T &val) {
  return "[" +
    ((val == numeric_limits<T>::max()) ? "infty" :
     ((val == -numeric_limits<T>::max()) ? "-infty" : to_string(val))) + "]";
}

template<class T> string
format_unsigned_like(T &val) {
  return "[" +
    ((val == numeric_limits<T>::max()) ? "infty" : to_string(val)) + "]";
}

template<class T> string
format_float_like(T &val) {
  return "[" +
    ((val == numeric_limits<T>::max()) ? "infty" :
     ((val == -numeric_limits<T>::max()) ? "-infty" :
      ((val == numeric_limits<T>::min()) ? "eps" :
       ((val == -numeric_limits<T>::min()) ? "-eps" :
        ((std::abs(val) < numeric_limits<T>::min()) ? "0.0" :
         to_string(val)))))) + "]";
}

string
Option::format_default_value() const {
  std::istringstream ss;
  if (arg_type == SMITHLAB_ARG_INT)
    return format_int_like(*int_value);
  else if (arg_type == SMITHLAB_ARG_LONG)
    return format_int_like(*long_value);
  else if (arg_type == SMITHLAB_ARG_UINT)
    return format_unsigned_like(*uint_value);
  else if (arg_type == SMITHLAB_ARG_ULONG)
    return format_unsigned_like(*ulong_value);
  else if (arg_type == SMITHLAB_ARG_FLOAT)
    return format_float_like(*float_value);
  else if (arg_type == SMITHLAB_ARG_DOUBLE)
    return format_float_like(*double_value);
  else if (arg_type == SMITHLAB_ARG_STRING)
    return *string_value;
  else if (arg_type == SMITHLAB_ARG_CHAR)
    return "[" + string(1, *char_value) + "]";
  else // if (arg_type == SMITHLAB_ARG_BOOL)
    return ""; //*bool_value ? "true" : "false";
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, int &val) :
  arg_type(SMITHLAB_ARG_INT), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), int_value(&val)  {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, unsigned int &val) :
  arg_type(SMITHLAB_ARG_UINT), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), uint_value(&val)  {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, long &val) :
  arg_type(SMITHLAB_ARG_LONG), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), long_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, unsigned long &val) :
  arg_type(SMITHLAB_ARG_ULONG), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), ulong_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, float &val) :
  arg_type(SMITHLAB_ARG_FLOAT), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), float_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, double &val) :
  arg_type(SMITHLAB_ARG_DOUBLE), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), double_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, string &val) :
  arg_type(SMITHLAB_ARG_STRING), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), string_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, bool &val) :
  arg_type(SMITHLAB_ARG_BOOL), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), bool_value(&val) {}

Option::Option(const string l_name, const char s_name, const string descr,
               const bool reqd, char &val) :
  arg_type(SMITHLAB_ARG_CHAR), long_name(l_name), short_name(s_name),
  description(descr), required(reqd), specified(false), char_value(&val) {}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

string
Option::format_option_name() const {
  std::ostringstream ss;
  if (short_name != '\0')
    ss << '-' << short_name << ", -" << long_name;
  else ss << "    -" << long_name;
  return ss.str();
}

string
Option::format_option_description(const size_t offset,
                                  const bool show_default) const {
  std::ostringstream ss;
  if (!description.empty()) {
    vector<string> parts;
    smithlab::split_whitespace(description, parts);
    if (required)
      parts.push_back("[required]");
    if (!required && show_default)
      parts.push_back(format_default_value());

    size_t line_len = 0;
    for (size_t i = 0; i < parts.size(); ++i) {
      if (offset + line_len + parts[i].size() >= MAX_LINE_LENGTH && i > 0) {
        line_len = 0;
        ss << endl;
      }
      if (i > 0 && line_len == 0)
        ss << string(offset, ' ');
      ss << parts[i] << " ";
      line_len += parts[i].size()+1; //+1 for the added space
    }
  }
  return ss.str();
}

bool
Option::option_match(const string &other) {
  return (long_name == other ||
          (other.length() > 1 && other[0] == '-' &&
           (other.substr(1) == long_name ||
            (other[1] == short_name && other.length() == 2))));
}

bool
Option::parse(vector<string> &command_line) {
  static const string dummy;
  if (!command_line.empty()) {
    for (size_t i = 0; i < command_line.size();)
      if (option_match(command_line[i])) {
        if (specified)
          throw runtime_error("duplicate use of option: " + long_name);

        if (i < command_line.size() - 1) {
          format_option(command_line[i + 1]);
        }
        else {
          // this will only work if it's a bool, because the
          // format_option function will ignore the argument
          format_option(dummy);
        }

        specified = true;
        // remove this option from the set of options
        command_line.erase(command_line.begin() + i);
        // if there was an argument (i.e. non bool) then remove that
        // argument also
        if (arg_type != SMITHLAB_ARG_BOOL) {
          command_line.erase(command_line.begin() + i);
        }
      }
      else {
        ++i;
      }
  }
  return (specified || !required);
}

void
Option::parse_config_file(vector<string> &options) {
  size_t i = 0;
  size_t op_num = options.size();
  while (i < op_num) {
    vector<string> opt_val = smithlab::split(options[i], "=");
    opt_val.front() = smithlab::strip(opt_val.front());
    opt_val.back() = smithlab::strip(opt_val.back());
    if (option_match(opt_val.front())) {
      format_option(opt_val.back());
      options.erase(options.begin() + i);
      specified = true;
      --op_num;
    }
    else {
      ++i;
    }
  }
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, int &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, unsigned &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, long &val)  {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, unsigned long &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, float &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, double &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, string &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, bool &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

void
OptionParser::add_opt(const string l_name, const char s_name, const string descr,
                      const bool reqd, char &val) {
  options.push_back(Option(l_name, s_name, descr, reqd, val));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

bool valid_option_char(char ch) {
  return std::isalnum(ch) || ch == '_';
}

static void
read_config_file(const string &config_filename,
                 vector<string> &config_file_options) {
  static const char comment_character = '#';
  static const char separator_character = ':';
  static const string outer_space = "^[:space:]+|[:space:]+$";
  static const string inner_space = "([:space:])[:space:]+";

  config_file_options.clear();

  std::ifstream in(config_filename);
  if (!in)
    throw runtime_error("cannot open config file: " + config_filename);

  string line;
  size_t line_number = 0;
  while (in) {

    if (!getline(in, line))
      throw runtime_error("failed to config line from " + config_filename);

    // remove leading and trailing space
    line = regex_replace(line, regex(outer_space), "");
    line = regex_replace(line, regex(inner_space), " ");

    if (!line.empty() && line.front() != comment_character) {

      const size_t sep_pos = line.find_first_of(separator_character);

      if (sep_pos == 0 || // catches ": "
          sep_pos >= line.length() - 1) // catches no sep or final char sep
        throw runtime_error("bad config file line: " + line);

      string option_label(line.substr(0, sep_pos));

      if (!all_of(begin(option_label), end(option_label), valid_option_char))
        throw runtime_error("bad option label: " + line);

      string option_value(line.substr(sep_pos + 1));
      // remove leading space
      option_value = regex_replace(line, regex(outer_space), "");

      if (!all_of(begin(option_value), end(option_value), valid_option_char))
        throw runtime_error("bad option label: " + line);

      // cerr << option_label << '\t' << option_value << endl;

      config_file_options.push_back(line);
    }
    in.peek();
    ++line_number;
  }
}

void
OptionParser::parse(const int argc, const char **argv,
                    vector<string> &arguments) {
  // The "2" below corresponds to the "about" and "help" options
  assert(options.size() >=  2);

  // The '1' and '+ 1' below is to skip over the program name
  arguments.clear();
  assert(argc >= 1);
  copy(argv + 1, argv + argc, back_inserter(arguments));

  // search for configuration file given in commnadline
  int i = 0;
  int arg_num = argc - 1;
  while (i < arg_num)
    if (arguments[i] == "--config") {
      vector<string> config_file_options;
      string config_filename;
      if (i + 1 < argc - 1)
        config_filename =  arguments[i+1];
      else
        // ads: need to check that this is really a filename
        throw runtime_error("--config requires config filename");
      read_config_file(config_filename, config_file_options);
      for (size_t j = 0; j < options.size(); ++j)
        options[j].parse_config_file(config_file_options);

      // ads: do we need to remove this arg? what if we need to know
      // that a config file was used?
      arguments.erase(arguments.begin() + i);
      arguments.erase(arguments.begin() + i);
      arg_num -= 2;
    }
    else
      ++i;

  // parse options given in commmand line
  for (size_t i = 0; i < options.size(); ++i)
    if (!options[i].parse(arguments) && first_missing_option_name.empty())
      first_missing_option_name = options[i].format_option_name();

  leftover_args = arguments;
}

void
OptionParser::parse(const int argc, const char **argv,
                    vector<string> &arguments, string config_filename) {
  // The "2" below corresponds to the "about" and "help" options
  assert(options.size() >=  2);

  if (!config_filename.empty()) {
    vector<string> config_file_options;
    read_config_file(config_filename, config_file_options);
    for (size_t i = 0; i < options.size(); ++i)
      options[i].parse_config_file(config_file_options);
  }
  arguments.clear();

  // The '1' and '+ 1' below is to skip over the program name
  assert(argc >= 1);
  copy(argv + 1, argv + argc, back_inserter(arguments));

  for (size_t i = 0; i < options.size(); ++i)
    if (!options[i].parse(arguments) && first_missing_option_name.empty())
      first_missing_option_name = options[i].format_option_name();

  leftover_args = arguments;
}

OptionParser::OptionParser(const string nm, const string descr,
                           string noflag_msg, const size_t n_left) :
  prog_name(nm), prog_descr(descr), noflag_message(noflag_msg),
  help_request(false), about_request(false),
  show_defaults(false), n_leftover(n_left) {
  add_opt("help", '?', "print this help message", false, help_request);
  add_opt("about", '\0', "print about message", false, about_request);
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//////
////// FOR PRINTING MESSAGES
//////

string
OptionParser::help_message() const {
  // corresponds to the two spaces before and
  static const string SPACE_BEFORE_SHORT = "  ";
  static const string SPACE_BTWN_SHRT_LNG = "  ";
  static const size_t TOTAL_ADDED_SPACE = 4;

  vector<string> option_names;
  size_t max_name_len = 0;
  for(size_t i = 0; i < options.size(); ++i) {
    option_names.push_back(options[i].format_option_name());
    max_name_len = std::max(max_name_len, option_names.back().length());
  }

  std::ostringstream ss;
  ss << "Usage: " << prog_name << " [OPTIONS]";
  if (!noflag_message.empty())
    ss << " " << noflag_message;
  ss << endl << endl;

  if (options.size() > 2) {
    ss << "Options:" << endl;
    // the loop below begins at 2 because the help and usage messages
    // are always the first two and treated separately
    for (size_t i = 2; i < options.size(); ++i)
      ss << SPACE_BEFORE_SHORT << std::left << std::setw(max_name_len)
         << option_names[i] << SPACE_BTWN_SHRT_LNG
         << options[i].format_option_description(max_name_len +
                                                 TOTAL_ADDED_SPACE,
                                                 show_defaults) << endl;
  }

  ss << endl << "Help options:" << endl;
  for (size_t i = 0; i < std::min(2ul, options.size()); ++i)
    ss << SPACE_BEFORE_SHORT << std::left << std::setw(max_name_len)
       << option_names[i] << SPACE_BTWN_SHRT_LNG
       << options[i].format_option_description(max_name_len +
                                               TOTAL_ADDED_SPACE,
                                               show_defaults) << endl;
  return ss.str();
}

string
OptionParser::about_message() const {
  static const char *PROGRAM_NAME_TAG = "PROGRAM: ";

  vector<string> parts;
  smithlab::split_whitespace(prog_descr, parts);

  std::ostringstream ss;
  ss << PROGRAM_NAME_TAG << prog_name << endl;
  ss << parts.front();
  size_t line_len = parts.front().length();
  for (size_t i = 1; i < parts.size(); ++i) {
    if (line_len + parts[i].size() >= MAX_LINE_LENGTH) {
      line_len = 0;
      ss << endl;
    }
    else ss << ' ';
    ss << parts[i];
    line_len += parts[i].length() + 1; // the "+1" is for the space
  }
  return ss.str();
}


string
OptionParser::invalid_leftover() const {
  static const string left_tag("invalid leftover args [should be ");
  static const string right_tag("]");

  std::ostringstream ss;
  if (n_leftover != std::numeric_limits<size_t>::max()) {
    ss << left_tag << n_leftover << right_tag << endl;
  }
  for (size_t i = 0; i < leftover_args.size(); ++i) {
    ss << "leftover arg #" << (i + 1) << "=\""
       << leftover_args[i] << "\"";
  }
  return ss.str();
}


string
OptionParser::option_missing_message() const {
  std::ostringstream ss;
  ss << "required argument missing: [" << first_missing_option_name << "]";
  return ss.str();
}
