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

#include "contaminants.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

// clang-format off
std::vector<std::pair<std::string, std::string>> contaminants = {  // NOLINT(cert-err58-cpp,cppcoreguidelines-avoid-non-const-global-variables)
  {"Illumina Single End Adapter 1", "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"},
  {"Illumina Single End Adapter 2", "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"},
  {"Illumina Single End PCR Primer 1", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Single End PCR Primer 2", "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"},
  {"Illumina Single End Sequencing Primer", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Paired End Adapter 1", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Paired End Adapter 2", "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"},
  {"Illumina Paried End PCR Primer 1", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Paired End PCR Primer 2", "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"},
  {"Illumina Paried End Sequencing Primer 1", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Paired End Sequencing Primer 2", "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"},
  {"Illumina DpnII expression Adapter 1", "ACAGGTTCAGAGTTCTACAGTCCGAC"},
  {"Illumina DpnII expression Adapter 2", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina DpnII expression PCR Primer 1", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina DpnII expression PCR Primer 2", "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
  {"Illumina DpnII expression Sequencing Primer", "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"},
  {"Illumina NlaIII expression Adapter 1", "ACAGGTTCAGAGTTCTACAGTCCGACATG"},
  {"Illumina NlaIII expression Adapter 2", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina NlaIII expression PCR Primer 1", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina NlaIII expression PCR Primer 2", "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
  {"Illumina NlaIII expression Sequencing Primer", "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"},
  {"Illumina Small RNA Adapter 1", "GTTCAGAGTTCTACAGTCCGACGATC"},
  {"Illumina Small RNA Adapter 2", "TGGAATTCTCGGGTGCCAAGG"},
  {"Illumina Small RNA RT Primer", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina Small RNA PCR Primer 2", "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
  {"Illumina Small RNA Sequencing Primer", "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"},
  {"Illumina Multiplexing Adapter 1", "GATCGGAAGAGCACACGTCT"},
  {"Illumina Multiplexing Adapter 2", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Multiplexing PCR Primer 1.01", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Multiplexing PCR Primer 2.01", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"},
  {"Illumina Multiplexing Read1 Sequencing Primer", "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"Illumina Multiplexing Index Sequencing Primer", "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"},
  {"Illumina Multiplexing Read2 Sequencing Primer", "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"},
  {"Illumina PCR Primer Index 1", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 2", "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 3", "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 4", "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 5", "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 6", "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 7", "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 8", "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 9", "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 10", "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 11", "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC"},
  {"Illumina PCR Primer Index 12", "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"},
  {"Illumina DpnII Gex Adapter 1", "GATCGTCGGACTGTAGAACTCTGAAC"},
  {"Illumina DpnII Gex Adapter 1.01", "ACAGGTTCAGAGTTCTACAGTCCGAC"},
  {"Illumina DpnII Gex Adapter 2", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina DpnII Gex Adapter 2.01", "TCGTATGCCGTCTTCTGCTTG"},
  {"Illumina DpnII Gex PCR Primer 1", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina DpnII Gex PCR Primer 2", "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
  {"Illumina DpnII Gex Sequencing Primer", "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"},
  {"Illumina NlaIII Gex Adapter 1.01", "TCGGACTGTAGAACTCTGAAC"},
  {"Illumina NlaIII Gex Adapter 1.02", "ACAGGTTCAGAGTTCTACAGTCCGACATG"},
  {"Illumina NlaIII Gex Adapter 2.01", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina NlaIII Gex Adapter 2.02", "TCGTATGCCGTCTTCTGCTTG"},
  {"Illumina NlaIII Gex PCR Primer 1", "CAAGCAGAAGACGGCATACGA"},
  {"Illumina NlaIII Gex PCR Primer 2", "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
  {"Illumina NlaIII Gex Sequencing Primer", "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"},
  {"Illumina 5p RNA Adapter", "GTTCAGAGTTCTACAGTCCGACGATC"},
  {"Illumina RNA Adapter1", "TGGAATTCTCGGGTGCCAAGG"},
  {"Illumina Small RNA 3p Adapter 1", "ATCTCGTATGCCGTCTTCTGCTTG"},
  {"Illumina Small RNA PCR Primer 1", "CAAGCAGAAGACGGCATACGA"},
  {"TruSeq Universal Adapter", "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
  {"TruSeq Adapter, Index 1", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 2", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 3", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 4", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 5", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 6", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 7", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 8", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 9", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 10", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 11", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 12", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 13", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 14", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 15", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 16", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 18", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 19", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 20", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 21", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 22", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 23", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 25", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG"},
  {"TruSeq Adapter, Index 27", "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG"},
  {"Illumina RNA RT Primer", "GCCTTGGCACCCGAGAATTCCA"},
  {"Illumina RNA PCR Primer", "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"},
  {"RNA PCR Primer, Index 1", "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 2", "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 3", "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 4", "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 5", "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 6", "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 7", "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 8", "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 9", "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 10", "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 11", "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 12", "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 13", "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 14", "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 15", "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 16", "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 17", "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 18", "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 19", "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 20", "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 21", "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 22", "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 23", "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 24", "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 25", "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 26", "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 27", "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 28", "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 29", "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 30", "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 31", "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 32", "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 33", "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 34", "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 35", "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 36", "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 37", "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 38", "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 39", "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 40", "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 41", "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 42", "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 43", "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 44", "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 45", "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 46", "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 47", "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"RNA PCR Primer, Index 48", "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
  {"ABI Dynabead EcoP Oligo", "CTGATCTAGAGGTACCGGATCCCAGCAGT"},
  {"ABI Solid3 Adapter A", "CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"},
  {"ABI Solid3 Adapter B", "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"},
  {"ABI Solid3 5' AMP Primer", "CCACTACGCCTCCGCTTTCCTCTCTATG"},
  {"ABI Solid3 3' AMP Primer", "CTGCCCCGGGTTCCTCATTCT"},
  {"ABI Solid3 EF1 alpha Sense Primer", "CATGTGTGTTGAGAGCTTC"},
  {"ABI Solid3 EF1 alpha Antisense Primer", "GAAAACCAAAGTGGTCCAC"},
  {"ABI Solid3 GAPDH Forward Primer", "TTAGCACCCCTGGCCAAGG"},
  {"ABI Solid3 GAPDH Reverse Primer", "CTTACTCCTTGGAGGCCATG"},
  {"Clontech Universal Primer Mix Short", "CTAATACGACTCACTATAGGGC"},
  {"Clontech Universal Primer Mix Long", "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"},
  {"Clontech SMARTer II A Oligonucleotide", "AAGCAGTGGTATCAACGCAGAGTAC"},
  {"Clontech SMART CDS Primer II A", "AAGCAGTGGTATCAACGCAGAGTACT"},
  {"Clontech_Universal_Primer_Mix_Short", "CTAATACGACTCACTATAGGGC"},
  {"Clontech_Universal_Primer_Mix_Long", "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"},
  {"Clontech_SMARTer_II_A_Oligonucleotide", "AAGCAGTGGTATCAACGCAGAGTAC"},
  {"Clontech_SMART_CDS_Primer_II_A", "AAGCAGTGGTATCAACGCAGAGTACT"},
  {"Clontech_SMART_CDS_Primer_II_A", "ACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGC"},
  {"Clontech_SMART_CDS_Primer_II_A", "GAGTACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGT"},
};
// clang-format on

auto
load_contaminants(const std::string &filename) -> void {
  // ADS: (todo) handle carriage returns and other control chars
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("failed to open contaminants file: " + filename);
  contaminants.clear();
  std::string line_data;
  while (std::getline(in, line_data)) {
    std::string_view line = line_data;
    const auto to_keep_prefix = line.find_first_not_of(" \t");
    if (to_keep_prefix == std::string_view::npos)
      continue;
    line.remove_prefix(std::min(to_keep_prefix, std::size(line)));
    if (line[0] == '#')
      continue;
    const auto to_keep_suffix = line.find_last_not_of(" \t");
    if (to_keep_suffix == std::string_view::npos)
      continue;
    line.remove_suffix(std::size(line) - to_keep_suffix - 1);
    std::string cleaned_line;
    for (auto itr = std::cbegin(line); itr != std::cend(line); ++itr)
      if (!std::isblank(*itr) ||
          (std::next(itr) != std::cend(line) && *itr != *std::next(itr)))
        cleaned_line += *itr;
    const auto tab_pos = cleaned_line.find('\t');
    if (tab_pos == std::string::npos ||
        tab_pos != cleaned_line.find_last_of('\t'))
      throw std::runtime_error("malformed line: " + line_data);
    const auto is_print = [](const auto c) { return std::isprint(c); };
    const auto name = cleaned_line.substr(0, tab_pos);
    const auto seq = cleaned_line.substr(tab_pos + 1);
    if (!std::ranges::all_of(name, is_print) ||
        !std::ranges::all_of(seq, is_print))
      throw std::runtime_error("malformed line: " + line_data);
    contaminants.emplace_back(name, seq);
  }
}

[[nodiscard]] auto
get_contam_name(const std::int64_t contam_idx) -> const std::string & {
  using std::string_literals::operator""s;
  static constexpr auto no_hit_label = "No Hit"s;
  if (contam_idx < 0 || contam_idx >= std::ssize(contaminants))
    return no_hit_label;
  return contaminants[contam_idx].first;
}

auto
load_contaminants(const std::string &filename) -> void;

// get the longest substring of left that is a prefix of right
[[nodiscard]] static inline auto
get_overlap(const auto &left, const auto &right) {
  const auto left_beg = std::cbegin(left);
  const auto left_end = std::cend(left);
  auto best_n_matches = 0L;
  for (auto left_itr = left_beg; left_itr != left_end; ++left_itr) {
    const auto [mm_left, _] =
      std::ranges::mismatch(std::ranges::subrange(left_itr, left_end), right);
    const auto n_matches = std::distance(left_itr, mm_left);
    best_n_matches = std::max(best_n_matches, n_matches);
  }
  return best_n_matches;
}

[[nodiscard]] auto
match_contaminant(const std::string &query) -> std::int64_t {
  auto best_idx = 0;
  auto best_match = 0L;
  auto best_match_len = 0LU;
  auto idx = 0;
  for (const auto &[name, seq] : contaminants) {
    const auto n_match =
      std::max(get_overlap(query, seq), get_overlap(seq, query));
    if (n_match > best_match) {
      best_idx = idx;
      best_match = n_match;
      best_match_len = std::size(seq);
    }
    ++idx;
  }
  const auto match_cutoff = std::min(best_match_len, std::size(query)) / 2.0;
  // If any sequence is a match, return the best one
  return best_match < match_cutoff ? -1 : best_idx;
}
