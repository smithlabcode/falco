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

#ifndef SRC_CONTAMINANTS_HPP_
#define SRC_CONTAMINANTS_HPP_

#include <array>
#include <string>
#include <utility>

// gets the largest suffix of left which is a prefix of right
// returns 0 if none exist
[[nodiscard]] static inline auto
get_overlap(const auto &left, const auto &right) {
  const auto right_beg = std::cbegin(right);
  const auto right_end = std::cend(right);
  const auto left_beg = std::cbegin(left);
  const auto left_end = std::cend(left);
  auto best_n_matches = 0l;
  for (auto left_itr = left_beg; left_itr != left_end; ++left_itr) {
    const auto mm_res = std::mismatch(left_itr, left_end, right_beg, right_end);
    // overlap must cover a suffix of left or the entirety of right
    if (mm_res.first == left_end || mm_res.second == right_end) {
      const auto n_matches = std::distance(left_itr, mm_res.first);
      best_n_matches = std::max(best_n_matches, n_matches);
    }
  }
  return best_n_matches;
}

[[nodiscard]] inline auto
match_contaminant(const auto &query, const auto &contams) -> std::string {
  static constexpr auto no_hit_label = "No Hit";
  std::string best_name;
  auto best_match = 0l;
  auto best_match_len = 0ul;
  for (const auto &[name, seq] : contams) {
    const auto n_match =
      std::max(get_overlap(query, seq), get_overlap(seq, query));
    if (n_match > best_match) {
      best_name = name;
      best_match = n_match;
      best_match_len = std::size(seq);
    }
  }
  const auto match_cutoff = std::min(best_match_len, std::size(query)) / 2.0;
  // If any sequence is a match, return the best one
  return best_match < match_cutoff ? no_hit_label : best_name;
}

using std::string_view_literals::operator""sv;

// clang-format off
static constexpr auto contaminants = std::array{
  std::pair{"Illumina Single End Adapter 1"sv, "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"Illumina Single End Adapter 2"sv, "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"sv},
  std::pair{"Illumina Single End PCR Primer 1"sv, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Single End PCR Primer 2"sv, "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"sv},
  std::pair{"Illumina Single End Sequencing Primer"sv, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Paired End Adapter 1"sv, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Paired End Adapter 2"sv, "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"sv},
  std::pair{"Illumina Paried End PCR Primer 1"sv, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Paired End PCR Primer 2"sv, "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"sv},
  std::pair{"Illumina Paried End Sequencing Primer 1"sv, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Paired End Sequencing Primer 2"sv, "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"sv},
  std::pair{"Illumina DpnII expression Adapter 1"sv, "ACAGGTTCAGAGTTCTACAGTCCGAC"sv},
  std::pair{"Illumina DpnII expression Adapter 2"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina DpnII expression PCR Primer 1"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina DpnII expression PCR Primer 2"sv, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"sv},
  std::pair{"Illumina DpnII expression Sequencing Primer"sv, "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"sv},
  std::pair{"Illumina NlaIII expression Adapter 1"sv, "ACAGGTTCAGAGTTCTACAGTCCGACATG"sv},
  std::pair{"Illumina NlaIII expression Adapter 2"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina NlaIII expression PCR Primer 1"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina NlaIII expression PCR Primer 2"sv, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"sv},
  std::pair{"Illumina NlaIII expression Sequencing Primer"sv, "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"sv},
  std::pair{"Illumina Small RNA Adapter 1"sv, "GTTCAGAGTTCTACAGTCCGACGATC"sv},
  std::pair{"Illumina Small RNA Adapter 2"sv, "TGGAATTCTCGGGTGCCAAGG"sv},
  std::pair{"Illumina Small RNA RT Primer"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina Small RNA PCR Primer 2"sv, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"sv},
  std::pair{"Illumina Small RNA Sequencing Primer"sv, "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"sv},
  std::pair{"Illumina Multiplexing Adapter 1"sv, "GATCGGAAGAGCACACGTCT"sv},
  std::pair{"Illumina Multiplexing Adapter 2"sv, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Multiplexing PCR Primer 1.01"sv, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Multiplexing PCR Primer 2.01"sv, "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"sv},
  std::pair{"Illumina Multiplexing Read1 Sequencing Primer"sv, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"Illumina Multiplexing Index Sequencing Primer"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"sv},
  std::pair{"Illumina Multiplexing Read2 Sequencing Primer"sv, "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"sv},
  std::pair{"Illumina PCR Primer Index 1"sv, "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 2"sv, "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 3"sv, "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 4"sv, "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 5"sv, "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 6"sv, "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 7"sv, "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 8"sv, "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 9"sv, "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 10"sv, "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 11"sv, "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC"sv},
  std::pair{"Illumina PCR Primer Index 12"sv, "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"sv},
  std::pair{"Illumina DpnII Gex Adapter 1"sv, "GATCGTCGGACTGTAGAACTCTGAAC"sv},
  std::pair{"Illumina DpnII Gex Adapter 1.01"sv, "ACAGGTTCAGAGTTCTACAGTCCGAC"sv},
  std::pair{"Illumina DpnII Gex Adapter 2"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina DpnII Gex Adapter 2.01"sv, "TCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"Illumina DpnII Gex PCR Primer 1"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina DpnII Gex PCR Primer 2"sv, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"sv},
  std::pair{"Illumina DpnII Gex Sequencing Primer"sv, "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"sv},
  std::pair{"Illumina NlaIII Gex Adapter 1.01"sv, "TCGGACTGTAGAACTCTGAAC"sv},
  std::pair{"Illumina NlaIII Gex Adapter 1.02"sv, "ACAGGTTCAGAGTTCTACAGTCCGACATG"sv},
  std::pair{"Illumina NlaIII Gex Adapter 2.01"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina NlaIII Gex Adapter 2.02"sv, "TCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"Illumina NlaIII Gex PCR Primer 1"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"Illumina NlaIII Gex PCR Primer 2"sv, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"sv},
  std::pair{"Illumina NlaIII Gex Sequencing Primer"sv, "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"sv},
  std::pair{"Illumina 5p RNA Adapter"sv, "GTTCAGAGTTCTACAGTCCGACGATC"sv},
  std::pair{"Illumina RNA Adapter1"sv, "TGGAATTCTCGGGTGCCAAGG"sv},
  std::pair{"Illumina Small RNA 3p Adapter 1"sv, "ATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"Illumina Small RNA PCR Primer 1"sv, "CAAGCAGAAGACGGCATACGA"sv},
  std::pair{"TruSeq Universal Adapter"sv, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"sv},
  std::pair{"TruSeq Adapter, Index 1"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 2"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 3"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 4"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 5"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 6"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 7"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 8"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 9"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 10"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 11"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 12"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 13"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 14"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 15"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 16"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 18"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 19"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 20"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 21"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 22"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 23"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 25"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"TruSeq Adapter, Index 27"sv, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG"sv},
  std::pair{"Illumina RNA RT Primer"sv, "GCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"Illumina RNA PCR Primer"sv, "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"sv},
  std::pair{"RNA PCR Primer, Index 1"sv, "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 2"sv, "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 3"sv, "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 4"sv, "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 5"sv, "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 6"sv, "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 7"sv, "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 8"sv, "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 9"sv, "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 10"sv, "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 11"sv, "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 12"sv, "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 13"sv, "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 14"sv, "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 15"sv, "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 16"sv, "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 17"sv, "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 18"sv, "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 19"sv, "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 20"sv, "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 21"sv, "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 22"sv, "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 23"sv, "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 24"sv, "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 25"sv, "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 26"sv, "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 27"sv, "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 28"sv, "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 29"sv, "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 30"sv, "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 31"sv, "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 32"sv, "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 33"sv, "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 34"sv, "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 35"sv, "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 36"sv, "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 37"sv, "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 38"sv, "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 39"sv, "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 40"sv, "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 41"sv, "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 42"sv, "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 43"sv, "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 44"sv, "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 45"sv, "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 46"sv, "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 47"sv, "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"RNA PCR Primer, Index 48"sv, "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"sv},
  std::pair{"ABI Dynabead EcoP Oligo"sv, "CTGATCTAGAGGTACCGGATCCCAGCAGT"sv},
  std::pair{"ABI Solid3 Adapter A"sv, "CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"sv},
  std::pair{"ABI Solid3 Adapter B"sv, "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"sv},
  std::pair{"ABI Solid3 5' AMP Primer"sv, "CCACTACGCCTCCGCTTTCCTCTCTATG"sv},
  std::pair{"ABI Solid3 3' AMP Primer"sv, "CTGCCCCGGGTTCCTCATTCT"sv},
  std::pair{"ABI Solid3 EF1 alpha Sense Primer"sv, "CATGTGTGTTGAGAGCTTC"sv},
  std::pair{"ABI Solid3 EF1 alpha Antisense Primer"sv, "GAAAACCAAAGTGGTCCAC"sv},
  std::pair{"ABI Solid3 GAPDH Forward Primer"sv, "TTAGCACCCCTGGCCAAGG"sv},
  std::pair{"ABI Solid3 GAPDH Reverse Primer"sv, "CTTACTCCTTGGAGGCCATG"sv},
  std::pair{"Clontech Universal Primer Mix Short"sv, "CTAATACGACTCACTATAGGGC"sv},
  std::pair{"Clontech Universal Primer Mix Long"sv, "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"sv},
  std::pair{"Clontech SMARTer II A Oligonucleotide"sv, "AAGCAGTGGTATCAACGCAGAGTAC"sv},
  std::pair{"Clontech SMART CDS Primer II A"sv, "AAGCAGTGGTATCAACGCAGAGTACT"sv},
  std::pair{"Clontech_Universal_Primer_Mix_Short"sv, "CTAATACGACTCACTATAGGGC"sv},
  std::pair{"Clontech_Universal_Primer_Mix_Long"sv, "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"sv},
  std::pair{"Clontech_SMARTer_II_A_Oligonucleotide"sv, "AAGCAGTGGTATCAACGCAGAGTAC"sv},
  std::pair{"Clontech_SMART_CDS_Primer_II_A"sv, "AAGCAGTGGTATCAACGCAGAGTACT"sv},
  std::pair{"Clontech_SMART_CDS_Primer_II_A"sv, "ACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGC"sv},
  std::pair{"Clontech_SMART_CDS_Primer_II_A"sv, "GAGTACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGT"sv},
};
// clang-format on

#endif  // SRC_CONTAMINANTS_HPP_
