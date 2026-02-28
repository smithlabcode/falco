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

#include <algorithm>
#include <chrono>
#include <fstream>
#include <sstream>

void
HtmlMaker::put_data(const std::string &placeholder, const std::string &data) {
  auto pos = html_boilerplate.find(placeholder);
  // Placeholder not found
  if (pos == std::string::npos)
    throw std::runtime_error("placeholder not found: " + placeholder);

  // at least one placeholder found
  while (pos != std::string::npos) {
    html_boilerplate.replace(pos, placeholder.size(), data);
    pos = html_boilerplate.find(placeholder, pos + 1);
  }
}

// Comments out html pieces if analyses were skipped
void
HtmlMaker::put_comment(std::string &comment_begin, std::string &comment_end,
                       bool done) {
  // put html comments if analysis was skipped
  if (!done) {
    put_data(comment_begin, "<!--");
    put_data(comment_end, "-->");
  }

  // otherwise delete placeholder
  else {
    put_data(comment_begin, "");
    put_data(comment_end, "");
  }
}

void
HtmlMaker::put_file_details(const FalcoConfig &falco_config) {
  // Put file name in filename placeholder
  put_data("{{filename}}", falco_config.filename_stripped);

  // Put date on date placeholder
  auto tmp =
    std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
  std::string time_fmt = std::string(ctime(&tmp));
  put_data("{{date}}", time_fmt);
}

HtmlMaker::HtmlMaker() { html_boilerplate = FalcoConfig::html_template; }
