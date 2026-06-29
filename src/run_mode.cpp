/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include "run_mode.hpp"

#include "boost/boost_unordered.hpp"

#include <iterator>
#include <string>
#include <vector>

// ADS: mapping between labels and function/var names

// clang-format off
//
// Label: what appears in config file
// var/fun name: how members of run_mode class are named
// Affect processing: if yes, then changes how analysis is done
//
// Label            | var/fun name  | affect processing
// ------------------------------------------------------------
// adapter          | do_adap       | Yes
// duplication      | do_dups       | Yes
// gc_sequence      | do_gc_content |
// kmer             | do_kmers      | Yes
// n_content        | do_n_content  |
// overrepresented  | do_overrep    | Yes (equivalent to dups)
// quality_base     | do_qual_base  |
// quality_sequence | do_qual_seq   |
// sequence         | do_sequence   |
// sequence_length  | do_length     |
// tile             | do_tiles      | Yes
// ------------------------------------------------------------
// clang-format on

// clang-format off
std::vector<std::string> run_mode::labels{  // NOLINT(cert-err58-cpp)
  "adapter",
  "duplication",
  "gc_sequence",
  "kmer",
  "n_content",
  "overrepresented",
  "quality_base",
  "quality_sequence",
  "sequence",
  "sequence_length",
  "tile",
};
// clang-format on

auto
run_mode::assign(const boost::unordered_flat_map<std::string, bool> &modes)
  -> void {
  static const auto set_mode = [&](const std::string &label, auto &the_mode) {
    const auto itr = modes.find(label);
    if (itr != std::cend(modes))
      the_mode = itr->second ? 1 : -1;
  };
  set_mode("adapter", do_adap_);
  set_mode("duplication", do_dups_);
  set_mode("gc_sequence", do_gc_content_);
  set_mode("kmer", do_kmers_);
  set_mode("sequence_length", do_length_);
  set_mode("n_content", do_n_content_);
  set_mode("overrepresented", do_overrep_);
  set_mode("quality_base", do_qual_base_);
  set_mode("quality_sequence", do_qual_seq_);
  set_mode("sequence", do_sequence_);
  set_mode("tile", do_tiles_);
}

auto
run_mode::set_unassigned() -> void {
  // clang-format off
  if (do_adap_ == 0) do_adap_ = do_adap_default;
  if (do_dups_ == 0) do_dups_ = do_dups_default;
  if (do_gc_content_ == 0) do_gc_content_ = do_gc_content_default;
  if (do_kmers_ == 0) do_kmers_ = do_kmers_default;
  if (do_n_content_ == 0) do_n_content_ = do_n_content_default;
  if (do_overrep_ == 0) do_overrep_ = do_overrep_default;
  if (do_qual_base_ == 0) do_qual_base_ = do_qual_base_default;
  if (do_qual_seq_ == 0) do_qual_seq_ = do_qual_seq_default;
  if (do_sequence_ == 0) do_sequence_ = do_sequence_default;
  if (do_length_ == 0) do_length_ = do_length_default;
  if (do_tiles_ == 0) do_tiles_ = do_tiles_default;
  // 'groups' not set in config file
  if (do_groups_ == 0) do_groups_ = do_groups_default;
  // clang-format on
}
