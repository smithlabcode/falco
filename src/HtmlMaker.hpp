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

#ifndef HTMLMAKER_HPP
#define HTMLMAKER_HPP

#include <string>
#include <fstream>

#include "FalcoConfig.hpp"
#include "FastqStats.hpp"

/*******************************************************/
/*************** HTML MAKER ****************************/
/*******************************************************/
class HtmlMaker {
public:
  explicit HtmlMaker(std::string html_template_path);
  // Fill data from module
  void put_data(const std::string &placeholder, const std::string &data);

  // Comment or remove placeholders
  void put_comment(std::string &comment_begin,
                   std::string &comment_end,
                   bool done);

  // Put file details and date
  void put_file_details(const FalcoConfig &falco_config);
  std::string html_boilerplate;
};
#endif
