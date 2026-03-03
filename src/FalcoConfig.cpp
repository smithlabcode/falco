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

#include <sys/stat.h>
#include <unistd.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

const std::string FalcoConfig::FalcoVersion = VERSION;

/*********************************************************/
/************** DEFAULT VALUES FOR FILES *****************/
/*********************************************************/

namespace FileConstants {
// These will become const bools in the stream reader
static const std::unordered_map<std::string,
                                std::unordered_map<std::string, double>>
  limits = {{"quality_base", {{"ignore", 0}}},
            {"duplication", {{"ignore", 0}, {"warn", 70}, {"error", 50}}},
            {"kmer", {{"ignore", 1}, {"warn", 2}, {"error", 5}}},
            {"n_content", {{"ignore", 0}, {"warn", 5}, {"error", 20}}},
            {"overrepresented", {{"ignore", 0}, {"warn", 0.1}, {"error", 1}}},
            {"quality_base_lower", {{"warn", 10}, {"error", 5}}},
            {"quality_base_median", {{"warn", 25}, {"error", 20}}},
            {"sequence", {{"ignore", 0}, {"warn", 10}, {"error", 20}}},
            {"gc_sequence", {{"ignore", 0}, {"warn", 15}, {"error", 30}}},
            {"quality_sequence", {{"ignore", 0}, {"warn", 27}, {"error", 20}}},
            {"tile", {{"ignore", 0}, {"warn", 5}, {"error", 10}}},
            {"sequence_length", {{"ignore", 0}, {"warn", 1}, {"error", 1}}},
            {"adapter", {{"ignore", 0}, {"warn", 5}, {"error", 10}}}};

/*************** CONTAMINANTS *****************/
// clang-format off
static const auto contaminants = std::vector<std::pair<std::string, std::string>>{
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

static const std::size_t adapter_size = 12;
};  // namespace FileConstants

/********************************************************************/
/**************************** AUX FUNCTIONS *************************/
/********************************************************************/

// Check if line is not a comment or newline
static inline bool
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

// Check if a std::string ends with another,
// to be use to figure out the file format
inline bool
endswith(const std::string &fn, const std::string &extn) {
  if (std::size(extn) > std::size(fn))
    return false;
  return std::equal(std::crbegin(extn), std::crend(extn), std::rbegin(fn));
}

/********************************************************************/
/******************** FALCOCONFIG FUNCTIONS *************************/
/********************************************************************/
// Default config properties
FalcoConfig::FalcoConfig(const int argc, char *argv[]) {
  casava = false;
  nanopore = false;
  nofilter = false;
  extract = false;
  nogroup = false;
  is_bisulfite = false;
  is_reverse_complement = false;
  read_step = 1;
  format = "";
  threads = 1;
  contaminants_file =
    std::string(PROGRAM_PATH) + "/Configuration/contaminant_list.txt";
  adapters_file = std::string(PROGRAM_PATH) + "/Configuration/adapter_list.txt";
  limits_file = std::string(PROGRAM_PATH) + "/Configuration/limits.txt";

  clean_zero_bytes(contaminants_file);
  clean_zero_bytes(adapters_file);
  clean_zero_bytes(limits_file);

  quiet = false;
  progress = false;
  tmpdir = ".";

  is_sam = false;
  is_bam = false;
  is_fastq = false;
  is_fastq_gz = false;

  std::ostringstream ost;
  for (int i = 0; i < argc; ++i) {
    if (i != 0)
      ost << " ";
    ost << std::string(argv[i]);
  }
  call = ost.str();
}

const std::vector<std::string> FalcoConfig::values_to_check = {
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
  "duplication",
  "kmer",
  "n_content",
  "overrepresented",
  "quality_base_lower",
  "quality_base_median",
  "sequence",
  "gc_sequence",
  "quality_sequence",
  "tile",
  "sequence_length",
  "adapter",
};

template <class T>
bool
check_if_not_ignored(const T &limits_map, const std::string &limit) {
  if (limits_map.find(limit) == end(limits_map))
    throw std::runtime_error("no instructions for limit " + limit);

  const auto the_limit = limits_map.find(limit)->second;
  if (the_limit.find("ignore") == end(the_limit))
    throw std::runtime_error("'ignore' option not set for limit " + limit);

  const bool ret = (the_limit.find("ignore")->second == 0.0);

  return ret;
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
  is_sam = false;
  is_bam = false;
  is_fastq_gz = false;
  is_fastq = false;

  if (format.empty()) {
    auto filename_lower(filename);
    to_lower(filename_lower);
    if (endswith(filename, ".sam") || endswith(filename, ".sam_mapped")) {
      is_sam = true;
    }
#ifdef USE_HTS
    else if (endswith(filename, ".bam") || endswith(filename, ".bam_mapped")) {
      is_bam = true;
    }
#endif
    else if (endswith(filename, ".fastq.gz") || endswith(filename, ".fq.gz")) {
      is_fastq_gz = true;
    }
    else if (endswith(filename, ".fastq") || endswith(filename, ".fq")) {
      is_fastq = true;
    }
  }
  else {
    if (format == "sam")
      is_sam = true;
#ifdef USE_HTS
    else if (format == "bam")
      is_bam = true;
#endif
    else if (format == "fq.gz" || format == "fastq.gz")
      is_fastq_gz = true;
    else if (format == "fq" || format == "fastq")
      is_fastq = true;
    else
      throw std::runtime_error("unrecognized file format: " + format);
  }
}

void
FalcoConfig::read_limits() {
  limits = FileConstants::limits;
  if (!std::filesystem::exists(limits_file)) {
    if (!quiet)
      std::cerr << "[limits]\tWARNING: using default limits because "
                << "limits file does not exist: " << limits_file << "\n";
  }
  else {
    std::ifstream in(limits_file);
    if (!in)
      throw std::runtime_error("problem opening limits file: " + limits_file);

    if (!quiet)
      std::cerr << "[limits]\tusing file " << limits_file << "\n";

    // Variables to parse lines
    std::string line, instruction;
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

      if (std::find(std::cbegin(values_to_check), std::cend(values_to_check),
                    limit) == std::cend(values_to_check))
        throw std::runtime_error("unknown limit option: " + limit);

      if (instruction != "warn" && instruction != "error" &&
          instruction != "ignore")
        throw std::runtime_error("unknown instruction for limit " + limit +
                                 ": " + instruction);

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

std::size_t
hash_adapter(const std::string &s) {
  std::size_t ans = 0;
  for (std::size_t i = 0; i < std::size(s); ++i) {
    if (s[i] != 'A' && s[i] != 'C' && s[i] != 'T' && s[i] != 'G')
      throw std::runtime_error("Bad adapter (non-ATGC characters): " + s);
    ans = (ans << 2) | actg_to_2bit(s[i]);
  }
  return ans;
}

void
FalcoConfig::read_adapters() {
  if (!std::filesystem::exists(adapters_file)) {
    if (!quiet)
      std::cerr << "[adapters]\tWARNING: using default adapters because "
                << "adapters file does not exist: " << adapters_file << "\n";

    adapter_names = FileConstants::adapter_names;
    adapter_seqs = FileConstants::adapter_seqs;

    adapter_hashes.clear();
    for (size_t i = 0; i < std::size(adapter_seqs); ++i)
      adapter_hashes.push_back(hash_adapter(adapter_seqs[i]));

    shortest_adapter_size = adapter_size = std::size(adapter_seqs[0]);
    return;
  }

  std::ifstream in(adapters_file);
  if (!in)
    throw std::runtime_error("problem opening adapters file: " + adapters_file);

  if (!quiet)
    std::cerr << "[adapters]\tusing file " << adapters_file << "\n";

  std::string line, token;
  std::vector<std::string> line_by_space;
  std::string adapter_name, adapter_seq;

  // The adapters file has a space separated name, and the last instance is
  // the biological sequence

  adapter_size = 0;
  adapter_names.clear();
  adapter_seqs.clear();
  adapter_hashes.clear();
  do_adapter_optimized = true;

  while (std::getline(in, line)) {
    if (is_content_line(line)) {
      if (std::size(adapter_names) > Constants::max_adapters) {
        throw std::runtime_error(
          "You are testing too many adapters. The maximum "
          "number is 128!");
      }
      adapter_name = "";
      adapter_seq = "";

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
        if (std::size(adapter_seq) < shortest_adapter_size) {
          shortest_adapter_size = std::size(adapter_seq);
        }
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
        << "contaminants file does not exist: " << contaminants_file << "\n";
    contaminants = FileConstants::contaminants;
    return;
  }
  std::ifstream in(contaminants_file);
  if (!in)
    throw std::runtime_error("problem opening contaminants file: " +
                             contaminants_file);

  if (!quiet)
    std::cerr << "[contaminants]\tusing file " << contaminants_file << "\n";
  std::vector<std::string> line_by_space;

  // The contaminants file has a space separated name, and the last
  // instance is the biological sequence
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
