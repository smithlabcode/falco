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

#include <string>
#include <vector>

// clang-format off
std::vector<std::string> run_mode::labels{
  "duplication",
  "kmer",
  "n_content",
  "overrepresented",
  "quality_base",
  "sequence",
  "gc_sequence",
  "quality_sequence",
  "tile",
  "sequence_length",
  "adapter",
};
// clang-format on

auto
run_mode::assign(const boost::unordered_flat_map<std::string, bool> &modes)
  -> void {
  static const auto set_mode = [&](const auto &label, auto &the_mode) {
    const auto itr = modes.find(label);
    if (itr != std::cend(modes))
      the_mode = itr->second;
  };
  set_mode("tile", do_tiles);
  set_mode("kmer", do_kmers);
  set_mode("duplication", do_dups);
}
