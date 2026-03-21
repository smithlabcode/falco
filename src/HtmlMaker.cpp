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

#include "HtmlMaker.hpp"

#include "config.h"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <regex>
#include <sstream>

void
HtmlMaker::put_data(const std::string &placeholder, const std::string &data) {
  auto pos = html_boilerplate.find(placeholder);
  // Placeholder not found
  if (pos == std::string::npos)
    throw std::runtime_error("placeholder not found: " + placeholder);

  // at least one placeholder found
  while (pos != std::string::npos) {
    html_boilerplate.replace(pos, std::size(placeholder), data);
    pos = html_boilerplate.find(placeholder, pos + 1);
  }
}

// Comments out html pieces if analyses were skipped
void
HtmlMaker::put_comment(std::string &comment_begin, std::string &comment_end,
                       const bool done) {
  if (!done) {  // put html comments if analysis was skipped
    put_data(comment_begin, "<!--");
    put_data(comment_end, "-->");
  }
  else {  // otherwise delete placeholder
    put_data(comment_begin, "");
    put_data(comment_end, "");
  }
}

void
HtmlMaker::put_file_details(const FalcoConfig &falco_config) {
  using namespace std::string_literals;
  static constexpr auto left_tag = "\\{\\{";
  static constexpr auto right_tag = "\\}\\}";

  const auto filename_formatted = falco_config.filename_stripped;
  std::regex filename_re(left_tag + "filename"s + right_tag);
  html_boilerplate =
    std::regex_replace(html_boilerplate, filename_re, filename_formatted);

  using system_clock = std::chrono::system_clock;
  auto time_unformatted = system_clock::to_time_t(system_clock::now());
  std::string time_formatted = std::string(ctime(&time_unformatted));

  std::regex date_re(left_tag + "date"s + right_tag);
  html_boilerplate =
    std::regex_replace(html_boilerplate, date_re, time_formatted);

  std::regex version_re(left_tag + "version"s + right_tag);
  html_boilerplate = std::regex_replace(html_boilerplate, version_re, VERSION);
}
