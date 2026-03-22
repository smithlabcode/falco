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

#include "FalcoConfig.hpp"
#include "FastqStats.hpp"
#include "html_template.hpp"

#include "config.h"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <unordered_set>

using std::string_literals::operator""s;

// Default values for files
namespace FileConstants {
// These will become const bools in the stream reader
// clang-format off
static const std::unordered_map<std::string,
                                std::unordered_map<std::string, double>>
  limits = {
  {"quality_base",
   {{"ignore", 0}}
  },
  {"duplication",
   {{"ignore", 0},
    {"warn", 70},
    {"error", 50}}
  },
  {"kmer",
   {{"ignore", 1},
    {"warn", 2},
    {"error", 5}}
  },
  {"n_content",
   {{"ignore", 0},
    {"warn", 5},
    {"error", 20}}
  },
  {"overrepresented",
   {{"ignore", 0},
    {"warn", 0.1},
    {"error", 1}}
  },
  {"quality_base_lower",
   {{"warn", 10},
    {"error", 5}}
  },
  {"quality_base_median",
   {{"warn", 25},
    {"error", 20}}
  },
  {"sequence",
   {{"ignore", 0},
    {"warn", 10},
    {"error", 20}}
  },
  {"gc_sequence",
   {{"ignore", 0},
    {"warn", 15},
    {"error", 30}}
  },
  {"quality_sequence",
   {{"ignore", 0},
    {"warn", 27},
    {"error", 20}}
  },
  {"tile",
   {{"ignore", 0},
    {"warn", 5},
    {"error", 10}}
  },
  {"sequence_length",
   {{"ignore", 0},
    {"warn", 1},
    {"error", 1}}
  },
  {"adapter",
   {{"ignore", 0},
    {"warn", 5},
    {"error", 10}}
  },
};
// clang-format on

/*************** CONTAMINANTS *****************/
// clang-format off
static const auto default_contaminants = {
  std::pair{"Illumina Single End Adapter 1"s, "GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"s},
  {"Illumina Single End Adapter 2"s, "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"s},
  {"Illumina Single End PCR Primer 1"s, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Single End PCR Primer 2"s, "CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"s},
  {"Illumina Single End Sequencing Primer"s, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Paired End Adapter 1"s, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Paired End Adapter 2"s, "GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"s},
  {"Illumina Paried End PCR Primer 1"s, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Paired End PCR Primer 2"s, "CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"s},
  {"Illumina Paried End Sequencing Primer 1"s, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Paired End Sequencing Primer 2"s, "CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"s},
  {"Illumina DpnII expression Adapter 1"s, "ACAGGTTCAGAGTTCTACAGTCCGAC"s},
  {"Illumina DpnII expression Adapter 2"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina DpnII expression PCR Primer 1"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina DpnII expression PCR Primer 2"s, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"s},
  {"Illumina DpnII expression Sequencing Primer"s, "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"s},
  {"Illumina NlaIII expression Adapter 1"s, "ACAGGTTCAGAGTTCTACAGTCCGACATG"s},
  {"Illumina NlaIII expression Adapter 2"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina NlaIII expression PCR Primer 1"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina NlaIII expression PCR Primer 2"s, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"s},
  {"Illumina NlaIII expression Sequencing Primer"s, "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"s},
  {"Illumina Small RNA Adapter 1"s, "GTTCAGAGTTCTACAGTCCGACGATC"s},
  {"Illumina Small RNA Adapter 2"s, "TGGAATTCTCGGGTGCCAAGG"s},
  {"Illumina Small RNA RT Primer"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina Small RNA PCR Primer 2"s, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"s},
  {"Illumina Small RNA Sequencing Primer"s, "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"s},
  {"Illumina Multiplexing Adapter 1"s, "GATCGGAAGAGCACACGTCT"s},
  {"Illumina Multiplexing Adapter 2"s, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Multiplexing PCR Primer 1.01"s, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Multiplexing PCR Primer 2.01"s, "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"s},
  {"Illumina Multiplexing Read1 Sequencing Primer"s, "ACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"Illumina Multiplexing Index Sequencing Primer"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"s},
  {"Illumina Multiplexing Read2 Sequencing Primer"s, "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"s},
  {"Illumina PCR Primer Index 1"s, "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 2"s, "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 3"s, "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 4"s, "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 5"s, "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 6"s, "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 7"s, "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 8"s, "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 9"s, "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 10"s, "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 11"s, "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC"s},
  {"Illumina PCR Primer Index 12"s, "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"s},
  {"Illumina DpnII Gex Adapter 1"s, "GATCGTCGGACTGTAGAACTCTGAAC"s},
  {"Illumina DpnII Gex Adapter 1.01"s, "ACAGGTTCAGAGTTCTACAGTCCGAC"s},
  {"Illumina DpnII Gex Adapter 2"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina DpnII Gex Adapter 2.01"s, "TCGTATGCCGTCTTCTGCTTG"s},
  {"Illumina DpnII Gex PCR Primer 1"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina DpnII Gex PCR Primer 2"s, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"s},
  {"Illumina DpnII Gex Sequencing Primer"s, "CGACAGGTTCAGAGTTCTACAGTCCGACGATC"s},
  {"Illumina NlaIII Gex Adapter 1.01"s, "TCGGACTGTAGAACTCTGAAC"s},
  {"Illumina NlaIII Gex Adapter 1.02"s, "ACAGGTTCAGAGTTCTACAGTCCGACATG"s},
  {"Illumina NlaIII Gex Adapter 2.01"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina NlaIII Gex Adapter 2.02"s, "TCGTATGCCGTCTTCTGCTTG"s},
  {"Illumina NlaIII Gex PCR Primer 1"s, "CAAGCAGAAGACGGCATACGA"s},
  {"Illumina NlaIII Gex PCR Primer 2"s, "AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"s},
  {"Illumina NlaIII Gex Sequencing Primer"s, "CCGACAGGTTCAGAGTTCTACAGTCCGACATG"s},
  {"Illumina 5p RNA Adapter"s, "GTTCAGAGTTCTACAGTCCGACGATC"s},
  {"Illumina RNA Adapter1"s, "TGGAATTCTCGGGTGCCAAGG"s},
  {"Illumina Small RNA 3p Adapter 1"s, "ATCTCGTATGCCGTCTTCTGCTTG"s},
  {"Illumina Small RNA PCR Primer 1"s, "CAAGCAGAAGACGGCATACGA"s},
  {"TruSeq Universal Adapter"s, "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"s},
  {"TruSeq Adapter, Index 1"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 2"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 3"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 4"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 5"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 6"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 7"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 8"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 9"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 10"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 11"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 12"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 13"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 14"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 15"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 16"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 18"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 19"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 20"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 21"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 22"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 23"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 25"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG"s},
  {"TruSeq Adapter, Index 27"s, "GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG"s},
  {"Illumina RNA RT Primer"s, "GCCTTGGCACCCGAGAATTCCA"s},
  {"Illumina RNA PCR Primer"s, "AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"s},
  {"RNA PCR Primer, Index 1"s, "CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 2"s, "CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 3"s, "CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 4"s, "CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 5"s, "CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 6"s, "CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 7"s, "CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 8"s, "CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 9"s, "CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 10"s, "CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 11"s, "CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 12"s, "CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 13"s, "CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 14"s, "CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 15"s, "CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 16"s, "CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 17"s, "CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 18"s, "CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 19"s, "CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 20"s, "CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 21"s, "CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 22"s, "CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 23"s, "CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 24"s, "CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 25"s, "CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 26"s, "CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 27"s, "CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 28"s, "CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 29"s, "CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 30"s, "CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 31"s, "CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 32"s, "CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 33"s, "CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 34"s, "CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 35"s, "CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 36"s, "CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 37"s, "CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 38"s, "CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 39"s, "CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 40"s, "CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 41"s, "CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 42"s, "CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 43"s, "CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 44"s, "CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 45"s, "CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 46"s, "CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 47"s, "CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"RNA PCR Primer, Index 48"s, "CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"s},
  {"ABI Dynabead EcoP Oligo"s, "CTGATCTAGAGGTACCGGATCCCAGCAGT"s},
  {"ABI Solid3 Adapter A"s, "CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"s},
  {"ABI Solid3 Adapter B"s, "CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"s},
  {"ABI Solid3 5' AMP Primer"s, "CCACTACGCCTCCGCTTTCCTCTCTATG"s},
  {"ABI Solid3 3' AMP Primer"s, "CTGCCCCGGGTTCCTCATTCT"s},
  {"ABI Solid3 EF1 alpha Sense Primer"s, "CATGTGTGTTGAGAGCTTC"s},
  {"ABI Solid3 EF1 alpha Antisense Primer"s, "GAAAACCAAAGTGGTCCAC"s},
  {"ABI Solid3 GAPDH Forward Primer"s, "TTAGCACCCCTGGCCAAGG"s},
  {"ABI Solid3 GAPDH Reverse Primer"s, "CTTACTCCTTGGAGGCCATG"s},
  {"Clontech Universal Primer Mix Short"s, "CTAATACGACTCACTATAGGGC"s},
  {"Clontech Universal Primer Mix Long"s, "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"s},
  {"Clontech SMARTer II A Oligonucleotide"s, "AAGCAGTGGTATCAACGCAGAGTAC"s},
  {"Clontech SMART CDS Primer II A"s, "AAGCAGTGGTATCAACGCAGAGTACT"s},
  {"Clontech_Universal_Primer_Mix_Short"s, "CTAATACGACTCACTATAGGGC"s},
  {"Clontech_Universal_Primer_Mix_Long"s, "CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"s},
  {"Clontech_SMARTer_II_A_Oligonucleotide"s, "AAGCAGTGGTATCAACGCAGAGTAC"s},
  {"Clontech_SMART_CDS_Primer_II_A"s, "AAGCAGTGGTATCAACGCAGAGTACT"s},
  {"Clontech_SMART_CDS_Primer_II_A"s, "ACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGC"s},
  {"Clontech_SMART_CDS_Primer_II_A"s, "GAGTACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGT"s},
};
// clang-format on

/*************** ADAPTERS *********************/
// Name (eg: Illumina Small RNA adapter)
static const auto adapter_names = std::vector<std::string>{
  "Illumina Universal Adapter",
  "Illumina Small RNA 3 prime Adapter",
  "Illumina Small RNA 5 prime Adapter",
  "Nextera Transposase Sequence",
  "SOLID Small RNA Adapter",
};

// Actual string sequence (eg: ATTGCCACA)
// clang-format off
static const auto adapter_seqs = std::vector<std::string>{
  "AGATCGGAAGAG",
  "TGGAATTCTCGG",
  "GATCGTCGGACT",
  "CTGTCTCTTATA",
  "CGCCTTGGCCGT",
};
// clang-format on

static constexpr std::size_t adapter_size = 12;
};  // namespace FileConstants

/********************************************************************/
/**************************** AUX FUNCTIONS *************************/
/********************************************************************/

// Check if line is not a comment or newline
[[nodiscard]] static inline bool
is_content_line(const std::string &line) {
  // ADS: not sure why this needs to require the line to have a length greater
  // than 1 instead of 0.
  return line.front() != '#' && std::size(line) > 1;
}

// This function is necessary for conda: zero bytes appear for whatever reason
// in the compile-time-resolved PROGRAM_PATH variable, and files are not read
// properly if these bytes are not removed
void
clean_zero_bytes(std::string &filename) {
  filename.erase(std::remove(std::begin(filename), std::end(filename), '\0'),
                 std::end(filename));
}

// Check if a string ends with another, to be use to figure out the file format
[[nodiscard]] static inline bool
endswith(const std::string &fn, const std::string &extn) {
  if (std::size(extn) > std::size(fn))
    return false;
  return std::equal(std::crbegin(extn), std::crend(extn), std::rbegin(fn));
}

// FalcoConfig functions

// Default config properties
FalcoConfig::FalcoConfig(const int argc, char *argv[]) {
  read_step = 1;
  threads = 1;
  contaminants_file =
    std::string(PROGRAM_PATH) + "/Configuration/contaminant_list.txt";
  adapters_file = std::string(PROGRAM_PATH) + "/Configuration/adapter_list.txt";
  limits_file = std::string(PROGRAM_PATH) + "/Configuration/limits.txt";

  clean_zero_bytes(contaminants_file);
  clean_zero_bytes(adapters_file);
  clean_zero_bytes(limits_file);

  tmpdir = ".";

  std::ostringstream ost;
  for (int i = 0; i < argc; ++i) {
    if (i != 0)
      ost << " ";
    ost << std::string(argv[i]);
  }
  call = ost.str();
}

const std::unordered_set<std::string> FalcoConfig::values_to_check{
  "duplication"s,
  "kmer"s,
  "n_content"s,
  "overrepresented"s,
  "quality_base"s,
  "sequence"s,
  "gc_sequence"s,
  "quality_sequence"s,
  "tile"s,
  "sequence_length"s,
  "adapter"s,
  "duplication"s,
  "kmer"s,
  "n_content"s,
  "overrepresented"s,
  "quality_base_lower"s,
  "quality_base_median"s,
  "sequence"s,
  "gc_sequence"s,
  "quality_sequence"s,
  "tile"s,
  "sequence_length"s,
  "adapter"s,
};

template <typename T>
[[nodiscard]] bool
check_if_not_ignored(const T &limits_map, const std::string &limit) {
  const auto itr_lim = limits_map.find(limit);
  if (itr_lim == std::cend(limits_map))
    throw std::runtime_error("no instructions for limit " + limit);
  const auto itr_ignore = itr_lim->second.find("ignore");
  if (itr_ignore == std::cend(itr_lim->second))
    throw std::runtime_error("'ignore' option not set for limit " + limit);
  return itr_ignore->second == 0.0;
}

void
FalcoConfig::setup() {
  // Now check for the file format (FASTQ/SAM/BAM, compressed or not)
  define_file_format();

  // Get filename without absolute path
  filename_stripped = std::filesystem::path(filename).filename().string();

  // read which modules to run and the cutoffs for pass/warn/fail
  read_limits();

  // Read files for appropriate modules
  if (do_adapter)
    read_adapters();
  if (do_overrepresented)
    read_contaminants_file();
}

void
FalcoConfig::define_file_format() {
  const auto to_lower = [](auto &s) {
    std::transform(std::cbegin(s), std::cend(s), std::begin(s),
                   [](const auto c) { return std::tolower(c); });
  };
  to_lower(format);

  // reset, important bececause the same FalcoConfig object is used across
  // possibly multiple input files

  if (format.empty()) {
    auto filename_lower(filename);
    to_lower(filename_lower);
    if (endswith(filename, ".sam") || endswith(filename, ".sam_mapped"))
      infile_type = infile_type_t::sam;
#ifdef USE_HTS
    else if (endswith(filename, ".bam") || endswith(filename, ".bam_mapped"))
      infile_type = infile_type_t::bam;
#endif
    else if (endswith(filename, ".fastq.gz") || endswith(filename, ".fq.gz"))
      infile_type = infile_type_t::fastq_gz;
    else if (endswith(filename, ".fastq") || endswith(filename, ".fq"))
      infile_type = infile_type_t::fastq;
    else
      throw std::runtime_error("unrecognized file format: " + filename);
  }
  else {
    if (format == "sam")
      infile_type = infile_type_t::sam;
#ifdef USE_HTS
    else if (format == "bam")
      infile_type = infile_type_t::bam;
#endif
    else if (format == "fq.gz" || format == "fastq.gz")
      infile_type = infile_type_t::fastq_gz;
    else if (format == "fq" || format == "fastq")
      infile_type = infile_type_t::fastq;
    else
      throw std::runtime_error("unrecognized file format: " + format);
  }
}

void
FalcoConfig::read_limits() {
  static const auto valid_instructions = std::unordered_set{
    "warn"s,
    "ignore"s,
    "error"s,
  };
  limits = FileConstants::limits;
  if (!std::filesystem::exists(limits_file)) {
    if (!quiet)
      std::cerr << "[limits]\tWARNING: using default limits because "
                << "limits file does not exist: " << limits_file << '\n';
  }
  else {
    std::ifstream in(limits_file);
    if (!in)
      throw std::runtime_error("problem opening limits file: " + limits_file);

    if (!quiet)
      std::cerr << "[limits]\tusing file " << limits_file << '\n';

    // Variables to parse lines
    std::string line;
    std::string instruction;
    double value{};
    while (std::getline(in, line)) {
      // Checks if the line has something to be parsed
      if (!is_content_line(line))
        continue;
      // Every line is a limit, warn/error/ignore and the value
      std::istringstream iss(line);
      std::string limit;
      if (!(iss >> limit >> instruction >> value))
        throw std::runtime_error("malformed limits line: \"" + line + "\"");

      if (values_to_check.count(limit) == 0)
        throw std::runtime_error("unknown limit option: " + limit);
      if (valid_instructions.count(instruction) == 0)
        throw std::runtime_error("unknown instruction for " + limit + ": " +
                                 instruction);
      limits[limit][instruction] = value;
    }
  }

  // Get data from config that tells us which analyses to skip
  do_duplication = check_if_not_ignored(limits, "duplication");
  do_kmer = check_if_not_ignored(limits, "kmer");
  do_n_content = check_if_not_ignored(limits, "n_content");
  do_overrepresented = check_if_not_ignored(limits, "overrepresented");
  do_quality_base = check_if_not_ignored(limits, "quality_base");
  do_sequence = check_if_not_ignored(limits, "sequence");
  do_gc_sequence = check_if_not_ignored(limits, "gc_sequence");
  do_quality_sequence = check_if_not_ignored(limits, "quality_sequence");
  do_tile = check_if_not_ignored(limits, "tile");
  do_adapter = check_if_not_ignored(limits, "adapter");
  do_sequence_length = check_if_not_ignored(limits, "sequence_length");
  do_adapter_optimized = false;
}

[[nodiscard]] static inline std::size_t
hash_adapter(const std::string &s) {
  std::size_t ans = 0;
  for (std::size_t i = 0; i < std::size(s); ++i) {
    const auto c = s[i];
    if (c != 'A' && c != 'C' && c != 'T' && c != 'G')
      throw std::runtime_error("Bad adapter (non-ATGC characters): " + s);
    ans = (ans << 2) | actg_to_2bit(c);
  }
  return ans;
}

void
FalcoConfig::read_adapters() {
  if (!std::filesystem::exists(adapters_file)) {
    if (!quiet)
      std::cerr << "[adapters]\tWARNING: using default adapters because "
                << "adapters file does not exist: " << adapters_file << '\n';
    adapter_names = FileConstants::adapter_names;
    adapter_seqs = FileConstants::adapter_seqs;
    adapter_hashes.clear();
    for (const auto &adap : adapter_seqs)
      adapter_hashes.push_back(hash_adapter(adap));
    shortest_adapter_size = adapter_size = std::size(adapter_seqs[0]);
    return;
  }

  std::ifstream in(adapters_file);
  if (!in)
    throw std::runtime_error("problem opening adapters file: " + adapters_file);

  if (!quiet)
    std::cerr << "[adapters]\tusing file " << adapters_file << '\n';

  std::string line;
  std::string token;
  std::vector<std::string> line_by_space;
  std::string adapter_name;
  std::string adapter_seq;

  // The adapters file has a space separated name, and the last instance is the
  // biological sequence

  adapter_size = 0;
  adapter_names.clear();
  adapter_seqs.clear();
  adapter_hashes.clear();
  do_adapter_optimized = true;

  while (std::getline(in, line)) {
    if (is_content_line(line)) {
      if (std::size(adapter_names) > Constants::max_adapters)
        throw std::runtime_error(
          "You are testing too many adapters. The maximum number is 128!");
      adapter_name.clear();
      adapter_seq.clear();

      line_by_space.clear();
      std::istringstream iss(line);
      while (iss >> token)
        line_by_space.push_back(token);

      if (std::size(line_by_space) > 1) {
        adapter_name = line_by_space[0];
        for (size_t i = 1; i < std::size(line_by_space) - 1; ++i)
          adapter_name += " " + line_by_space[i];

        adapter_seq = line_by_space.back();
        if (std::size(adapter_seq) > 32) {
          std::cerr << "[adapters]\tadapter size is more then 32. Use slow "
                       "adapters search\n";
          do_adapter_optimized = false;
        }
      }

      // store information
      adapter_names.push_back(adapter_name);
      adapter_seqs.push_back(adapter_seq);
      adapter_hashes.push_back(hash_adapter(adapter_seq));

      if (adapter_size == 0) {
        adapter_size = std::size(adapter_seq);
        shortest_adapter_size = adapter_size;
      }
      else if (std::size(adapter_seq) != adapter_size) {
        std::cerr << "[adapters]\tadapters have different size. Use slow "
                     "adapters search\n";
        do_adapter_optimized = false;
        if (std::size(adapter_seq) < shortest_adapter_size)
          shortest_adapter_size = std::size(adapter_seq);
      }
    }
  }
}

void
FalcoConfig::read_contaminants_file() {
  if (!std::filesystem::exists(contaminants_file)) {
    if (!quiet)
      std::cerr
        << "[contaminants]\tWARNING: using default contaminants because "
        << "contaminants file does not exist: " << contaminants_file << '\n';
    contaminants = FileConstants::default_contaminants;
    return;
  }
  std::ifstream in(contaminants_file);
  if (!in)
    throw std::runtime_error("problem opening contaminants file: " +
                             contaminants_file);
  if (!quiet)
    std::cerr << "[contaminants]\tusing file " << contaminants_file << '\n';

  // The contaminants file has a space separated name, and the last instance is
  // the biological sequence
  std::vector<std::string> line_by_space;
  std::string line;
  contaminants.clear();
  while (std::getline(in, line)) {
    if (is_content_line(line)) {
      std::istringstream iss(line);
      std::string token;
      while (iss >> token)
        line_by_space.push_back(token);
      if (std::size(line_by_space) > 1) {
        std::string contaminant_name = line_by_space[0];
        for (size_t i = 1; i < std::size(line_by_space) - 1; ++i)
          contaminant_name += " " + line_by_space[i];
        contaminants.push_back(
          std::make_pair(contaminant_name, line_by_space.back()));
      }
      line_by_space.clear();
    }
  }
}

const std::string FalcoConfig::html_template = falco_config_html_template;
