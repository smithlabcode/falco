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
  explicit HtmlMaker(std::string filepath);
  void write(FastqStats &stats, FalcoConfig &fc);
  std::string sourcecode;

private:
  void replace_placeholder_with_data(const std::string &placeholder,
                                     const std::string &data);

  // Replace placeholders with html comments if analysis was skipped
  void comment_if_skipped(std::string ph_begin, std::string ph_end, bool skip);

  // First thing to do: file details
  void put_file_details(const FalcoConfig &fc);

  // Comment out analyses that were skipped
  void comment_out(const FalcoConfig &fc);

  // Functions to replace pass warn fail
  void put_pass_warn_fail(const FastqStats &stats);

  // Function to replace template placeholders with data
  void make_basic_statistics(const FastqStats &stats,
                             FalcoConfig &fc);

  void make_position_quality_data(const FastqStats &stats,
                                  const FalcoConfig &fc);

  void make_tile_quality_data(FastqStats &stats,
                              const FalcoConfig &fc);

  void make_sequence_quality_data(const FastqStats &stats,
                                  const FalcoConfig &fc);
  void make_base_sequence_content_data(const FastqStats &stats,
                                       const FalcoConfig &fc);

  void make_sequence_gc_content_data(const FastqStats &stats,
                                     const FalcoConfig &fc);

  void make_base_n_content_data(const FastqStats &stats,
                                const FalcoConfig &fc);

  void make_sequence_length_data(const FastqStats &stats,
                                 const FalcoConfig &fc);

  void make_sequence_duplication_data(const FastqStats &stats,
                                      const FalcoConfig &fc);

  void make_overrepresented_sequences_data(const FastqStats &stats,
                                           const FalcoConfig &fc);

  void make_adapter_content_data(FastqStats &stats,
                                 FalcoConfig &fc);
};
#endif
