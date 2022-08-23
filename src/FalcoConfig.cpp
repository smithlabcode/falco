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

#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <algorithm>
#include <sstream>

using std::ostringstream;
using std::transform;
using std::string;
using std::vector;
using std::unordered_map;
using std::pair;
using std::make_pair;
using std::ifstream;
using std::runtime_error;
using std::istringstream;
using std::cerr;

const string FalcoConfig::FalcoVersion = "1.0.0";

/*********************************************************/
/************** DEFAULT VALUES FOR FILES *****************/
/*********************************************************/

namespace FileConstants {
  // These will become const bools in the stream reader
  static const std::unordered_map<std::string,
               std::unordered_map<std::string, double> >
  limits = {
    {"quality_base",{{"ignore",0}}},
    {"duplication",{{"ignore",0}, {"warn",70}, {"error",50}}},
    {"kmer",{{"ignore", 1}, {"warn",2}, {"error",5}}},
    {"n_content",{{"ignore", 0}, {"warn",5}, {"error",20}}},
    {"overrepresented",{{"ignore", 0}, {"warn",0.1}, {"error",1}}},
    {"quality_base_lower",{{"warn",10}, {"error",5}}},
    {"quality_base_median",{{"warn",25}, {"error",20}}},
    {"sequence",{{"ignore", 0}, {"warn",10}, {"error",20}}},
    {"gc_sequence",{{"ignore", 0}, {"warn",15}, {"error",30}}},
    {"quality_sequence",{{"ignore", 0}, {"warn",27}, {"error",20}}},
    {"tile",{{"ignore", 0}, {"warn",5}, {"error",10}}},
    {"sequence_length",{{"ignore", 0}, {"warn",1}, {"error",1}}},
    {"adapter",{{"ignore", 0}, {"warn",5}, {"error",10}}}
  };

  /*************** CONTAMINANTS *****************/
  static const std::vector<std::pair<std::string, std::string> >
  contaminants = {
    {"Illumina Single End Adapter 1","GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"},
    {"Illumina Single End Adapter 2","CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"},
    {"Illumina Single End PCR Primer 1","AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Single End PCR Primer 2","CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT"},
    {"Illumina Single End Sequencing Primer","ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Paired End Adapter 1","ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Paired End Adapter 2","GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG"},
    {"Illumina Paried End PCR Primer 1","AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Paired End PCR Primer 2","CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"},
    {"Illumina Paried End Sequencing Primer 1","ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Paired End Sequencing Primer 2","CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT"},
    {"Illumina DpnII expression Adapter 1","ACAGGTTCAGAGTTCTACAGTCCGAC"},
    {"Illumina DpnII expression Adapter 2","CAAGCAGAAGACGGCATACGA"},
    {"Illumina DpnII expression PCR Primer 1","CAAGCAGAAGACGGCATACGA"},
    {"Illumina DpnII expression PCR Primer 2","AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
    {"Illumina DpnII expression Sequencing Primer","CGACAGGTTCAGAGTTCTACAGTCCGACGATC"},
    {"Illumina NlaIII expression Adapter 1","ACAGGTTCAGAGTTCTACAGTCCGACATG"},
    {"Illumina NlaIII expression Adapter 2","CAAGCAGAAGACGGCATACGA"},
    {"Illumina NlaIII expression PCR Primer 1","CAAGCAGAAGACGGCATACGA"},
    {"Illumina NlaIII expression PCR Primer 2","AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
    {"Illumina NlaIII expression Sequencing Primer","CCGACAGGTTCAGAGTTCTACAGTCCGACATG"},
    {"Illumina Small RNA Adapter 1","GTTCAGAGTTCTACAGTCCGACGATC"},
    {"Illumina Small RNA Adapter 2","TGGAATTCTCGGGTGCCAAGG"},
    {"Illumina Small RNA RT Primer","CAAGCAGAAGACGGCATACGA"},
    {"Illumina Small RNA PCR Primer 2","AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
    {"Illumina Small RNA Sequencing Primer","CGACAGGTTCAGAGTTCTACAGTCCGACGATC"},
    {"Illumina Multiplexing Adapter 1","GATCGGAAGAGCACACGTCT"},
    {"Illumina Multiplexing Adapter 2","ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Multiplexing PCR Primer 1.01","AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Multiplexing PCR Primer 2.01","GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"},
    {"Illumina Multiplexing Read1 Sequencing Primer","ACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"Illumina Multiplexing Index Sequencing Primer","GATCGGAAGAGCACACGTCTGAACTCCAGTCAC"},
    {"Illumina Multiplexing Read2 Sequencing Primer","GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT"},
    {"Illumina PCR Primer Index 1","CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 2","CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 3","CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 4","CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 5","CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 6","CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 7","CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 8","CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 9","CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 10","CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 11","CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC"},
    {"Illumina PCR Primer Index 12","CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC"},
    {"Illumina DpnII Gex Adapter 1","GATCGTCGGACTGTAGAACTCTGAAC"},
    {"Illumina DpnII Gex Adapter 1.01","ACAGGTTCAGAGTTCTACAGTCCGAC"},
    {"Illumina DpnII Gex Adapter 2","CAAGCAGAAGACGGCATACGA"},
    {"Illumina DpnII Gex Adapter 2.01","TCGTATGCCGTCTTCTGCTTG"},
    {"Illumina DpnII Gex PCR Primer 1","CAAGCAGAAGACGGCATACGA"},
    {"Illumina DpnII Gex PCR Primer 2","AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
    {"Illumina DpnII Gex Sequencing Primer","CGACAGGTTCAGAGTTCTACAGTCCGACGATC"},
    {"Illumina NlaIII Gex Adapter 1.01","TCGGACTGTAGAACTCTGAAC"},
    {"Illumina NlaIII Gex Adapter 1.02","ACAGGTTCAGAGTTCTACAGTCCGACATG"},
    {"Illumina NlaIII Gex Adapter 2.01","CAAGCAGAAGACGGCATACGA"},
    {"Illumina NlaIII Gex Adapter 2.02","TCGTATGCCGTCTTCTGCTTG"},
    {"Illumina NlaIII Gex PCR Primer 1","CAAGCAGAAGACGGCATACGA"},
    {"Illumina NlaIII Gex PCR Primer 2","AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA"},
    {"Illumina NlaIII Gex Sequencing Primer","CCGACAGGTTCAGAGTTCTACAGTCCGACATG"},
    {"Illumina 5p RNA Adapter","GTTCAGAGTTCTACAGTCCGACGATC"},
    {"Illumina RNA Adapter1","TGGAATTCTCGGGTGCCAAGG"},
    {"Illumina Small RNA 3p Adapter 1","ATCTCGTATGCCGTCTTCTGCTTG"},
    {"Illumina Small RNA PCR Primer 1","CAAGCAGAAGACGGCATACGA"},
    {"TruSeq Universal Adapter","AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT"},
    {"TruSeq Adapter, Index 1","GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 2","GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 3","GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 4","GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 5","GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 6","GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 7","GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 8","GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 9","GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 10","GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 11","GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 12","GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 13","GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 14","GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 15","GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 16","GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 18","GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 19","GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 20","GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 21","GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 22","GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 23","GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 25","GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG"},
    {"TruSeq Adapter, Index 27","GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG"},
    {"Illumina RNA RT Primer","GCCTTGGCACCCGAGAATTCCA"},
    {"Illumina RNA PCR Primer","AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA"},
    {"RNA PCR Primer, Index 1","CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 2","CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 3","CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 4","CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 5","CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 6","CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 7","CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 8","CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 9","CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 10","CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 11","CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 12","CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 13","CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 14","CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 15","CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 16","CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 17","CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 18","CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 19","CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 20","CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 21","CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 22","CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 23","CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 24","CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 25","CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 26","CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 27","CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 28","CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 29","CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 30","CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 31","CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 32","CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 33","CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 34","CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 35","CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 36","CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 37","CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 38","CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 39","CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 40","CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 41","CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 42","CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 43","CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 44","CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 45","CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 46","CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 47","CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"RNA PCR Primer, Index 48","CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"},
    {"ABI Dynabead EcoP Oligo","CTGATCTAGAGGTACCGGATCCCAGCAGT"},
    {"ABI Solid3 Adapter A","CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG"},
    {"ABI Solid3 Adapter B","CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT"},
    {"ABI Solid3 5' AMP Primer","CCACTACGCCTCCGCTTTCCTCTCTATG"},
    {"ABI Solid3 3' AMP Primer","CTGCCCCGGGTTCCTCATTCT"},
    {"ABI Solid3 EF1 alpha Sense Primer","CATGTGTGTTGAGAGCTTC"},
    {"ABI Solid3 EF1 alpha Antisense Primer","GAAAACCAAAGTGGTCCAC"},
    {"ABI Solid3 GAPDH Forward Primer","TTAGCACCCCTGGCCAAGG"},
    {"ABI Solid3 GAPDH Reverse Primer","CTTACTCCTTGGAGGCCATG"},
    {"Clontech Universal Primer Mix Short","CTAATACGACTCACTATAGGGC"},
    {"Clontech Universal Primer Mix Long","CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"},
    {"Clontech SMARTer II A Oligonucleotide","AAGCAGTGGTATCAACGCAGAGTAC"},
    {"Clontech SMART CDS Primer II A","AAGCAGTGGTATCAACGCAGAGTACT"},
    {"Clontech_Universal_Primer_Mix_Short","CTAATACGACTCACTATAGGGC"},
    {"Clontech_Universal_Primer_Mix_Long","CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT"},
    {"Clontech_SMARTer_II_A_Oligonucleotide","AAGCAGTGGTATCAACGCAGAGTAC"},
    {"Clontech_SMART_CDS_Primer_II_A","AAGCAGTGGTATCAACGCAGAGTACT"},
    {"Clontech_SMART_CDS_Primer_II_A","ACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGC"},
    {"Clontech_SMART_CDS_Primer_II_A","GAGTACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGT"}
  };

  /*************** ADAPTERS *********************/
  // Name (eg: Illumina Small RNA adapter)
  static const std::vector<std::string>
  adapter_names = {
    "Illumina Universal Adapter",
    "Illumina Small RNA 3 prime Adapter",
    "Illumina Small RNA 5 prime Adapter",
    "Nextera Transposase Sequence",
    "SOLID Small RNA Adapter",
  };

  // Actual string sequence (eg: ATTGCCACA)
  static const std::vector<std::string>
  adapter_seqs = {
    "AGATCGGAAGAG",
    "TGGAATTCTCGG",
    "GATCGTCGGACT",
    "CTGTCTCTTATA",
    "CGCCTTGGCCGT",
  };

  static const size_t adapter_size = 12;
};


/********************************************************************/
/**************************** AUX FUNCTIONS *************************/
/********************************************************************/

// Check if line is not a comment or newline
inline bool
is_content_line (const string &line) {
  // comment
  if (line[0] == '#')
    return false;

  // newline
  if (line.size() <= 1)
    return false;

  return true;
}

// Check existance of config files
inline bool
file_exists(const std::string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}

// Check if a std::string ends with another,
// to be use to figure out the file format
inline bool
endswith(std::string const &value, std::string const &ending) {
  if (ending.size() > value.size()) {
    return false;
  }
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

// Removes absolute path from a file
static string
strip_path(string full_path) {
  size_t start = full_path.find_last_of('/');
  if (start == string::npos)
    start = 0;
  else
    ++start;
  return full_path.substr(start);
}

/********************************************************************/
/******************** FALCOCONFIG FUNCTIONS *************************/
/********************************************************************/
// Default config properties
FalcoConfig::FalcoConfig(const int argc, const char **argv) {
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
  contaminants_file = string(PROGRAM_PATH) + "/Configuration/contaminant_list.txt";
  adapters_file = string(PROGRAM_PATH) + "/Configuration/adapter_list.txt";
  limits_file = string(PROGRAM_PATH) + "/Configuration/limits.txt";

  if (!file_exists(adapters_file))
    adapters_file = "";

  if (!file_exists(contaminants_file))
    contaminants_file = "";

  if (!file_exists(limits_file))
    limits_file = "";

  quiet = false;
  tmpdir = ".";

  is_sam = false;
  is_bam = false;
  is_fastq = false;
  is_fastq_gz = false;

  ostringstream ost;
  for (int i = 0; i < argc; ++i) {
    if (i != 0)
      ost << " " ;
    ost << string(argv[i]);
  }
  call = ost.str();
}

const vector<string> FalcoConfig::values_to_check({
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
    "adapter"
  });

void
FalcoConfig::setup() {
  // Now check for the file format (FASTQ/SAM/BAM, compressed or not)
  define_file_format();

  // Get filename without absolute path
  filename_stripped = strip_path(filename);

  // read which modules to run and the cutoffs for pass/warn/fail
  read_limits();
  // Read files for appropriate modules
  if (limits["adapter"]["ignore"] == 0.0)
    read_adapters();
  if (limits["adapter"]["ignore"] == 0.0)
    read_contaminants_file();
}

void
FalcoConfig::define_file_format() {
  transform(begin(format), end(format), begin(format), tolower);
  string tmp_filename = filename;
  transform(begin(tmp_filename), end(tmp_filename), begin(tmp_filename), tolower);
  if (format == "") {
    if (endswith(tmp_filename, "sam")) {
      is_sam = true;
    }
#ifdef USE_HTS
    if (endswith(tmp_filename, "bam")) {
      is_bam = true;
    }
#endif
    if (endswith(tmp_filename, "fastq.gz")) {
      is_fastq_gz = true;
    }
    if (endswith(tmp_filename, "fq.gz")) {
      is_fastq_gz = true;
    }
    if (endswith(tmp_filename, "fastq")) {
      is_fastq = true;
    }
    if (endswith(tmp_filename, "fq")) {
      is_fastq = true;
    }
  }
  else {
    if (format == "sam") is_sam = true;
#ifdef USE_HTS
    else if (format == "bam") is_bam = true;
#endif
    else if (format == "fq.gz" || format == "fastq.gz") is_fastq_gz = true;
    else if (format == "fq" || format == "fastq") is_fastq_gz = true;
    else throw runtime_error("unrecognized file format: " + format);
  }
}

void
FalcoConfig::read_limits() {
  if (limits_file == "") {
    if (!quiet)
      cerr << "[limits]\tusing default limit cutoffs (no file specified)\n";
    limits = FileConstants::limits;
  } else {
    ifstream in(limits_file);
    if (!in)
      throw runtime_error("limits file does not exist: " + limits_file);

    if (!quiet)
      cerr << "[limitst]\tusing file " << limits_file << "\n"; 

    // Variables to parse lines
    string line, instruction;
    double value;
    while (getline(in, line)) {
      // Checks if the line has something to be parsed
      if (is_content_line (line)) {
        istringstream iss(line);

        // Every line is a limit, warn/error/ignore and the value
        string limit;
        if (!(iss >> limit >> instruction >> value))
          throw runtime_error("malformed limits line: \"" + line + "\"");

        if (find(begin(values_to_check), end(values_to_check), limit)
            == end(values_to_check))
          throw runtime_error("unknown limit option: " + limit);

        if (instruction != "warn" &&
            instruction != "error" &&
            instruction != "ignore")
          throw runtime_error("unknown instruction for limit " +
                              limit + ": " + instruction);

        limits[limit][instruction] = value;
      }
    }
  }
  for (auto v : values_to_check)
    if (limits.count(v) == 0)
      throw runtime_error("instruction for limit " + v +
                          " not found in file " + limits_file);
  // Get useful data from config that tells us which analyses to skip
  do_duplication = (limits["duplication"]["ignore"] == 0.0);
  do_kmer = (limits["kmer"]["ignore"] == 0.0);
  do_n_content = (limits["n_content"]["ignore"] == 0.0);
  do_overrepresented = (limits["overrepresented"]["ignore"] == 0.0);
  do_quality_base = (limits["quality_base"]["ignore"] == 0.0);
  do_sequence = (limits["sequence"]["ignore"] == 0.0);
  do_gc_sequence = (limits["gc_sequence"]["ignore"] == 0.0);
  do_quality_sequence= (limits["quality_sequence"]["ignore"] == 0.0);
  do_tile = (limits["tile"]["ignore"] == 0.0);
  do_adapter = (limits["adapter"]["ignore"] == 0.0);
  do_sequence_length = (limits["sequence_length"]["ignore"] == 0.0);
}

size_t
hash_adapter(const string &s) {
  size_t ans = 0;
  for (size_t i = 0; i < s.size(); ++i) {
    if (s[i] != 'A' && s[i] != 'C' && s[i] != 'T' && s[i] != 'G')
      throw runtime_error("Bad adapter (non-ATGC characters): " + s);

    ans = (ans << 2) | actg_to_2bit(s[i]);
  }

  return ans;
}

void
FalcoConfig::read_adapters() {
  if (adapters_file == "") {
    if (!quiet)
      cerr << "[adapters]\tusing default adapters (no file specified)\n";
    adapter_names = FileConstants::adapter_names;
    adapter_seqs = FileConstants::adapter_seqs;

    adapter_hashes.clear();
    for (size_t i = 0; i < adapter_seqs.size(); ++i)
      adapter_hashes.push_back(hash_adapter(adapter_seqs[i]));

    adapter_size = adapter_seqs[0].size();
    return;
  }
  ifstream in(adapters_file);
  if (!in)
    throw runtime_error("adapter file not found: " + adapters_file);

  if (!quiet)
    cerr << "[adapters]\tusing file " << adapters_file << "\n";
  string line, _tmp;
  vector<string> line_by_space;
  string adapter_name, adapter_seq;

  // The adapters file has a space separated name, and the last instance is
  // the biological sequence

  adapter_size = 0;
  adapter_names.clear();
  adapter_seqs.clear();
  adapter_hashes.clear();
  do_adapter_optimized = true;

  while (getline(in, line)) {
    if (is_content_line(line)) {
      if (adapter_names.size() > Constants::max_adapters)
        throw runtime_error("You are testing too many adapters. The maximum "
                            "number is 128!");
      adapter_name = "";
      adapter_seq = "";
      istringstream iss(line);
      while (iss >> _tmp) {
        line_by_space.push_back(_tmp);
      }

      if (line_by_space.size() > 1) {
        for (size_t i = 0; i < line_by_space.size() - 1; ++i)
          adapter_name += line_by_space[i] + " ";
        adapter_seq = line_by_space.back();

        if (adapter_seq.size() > 32) {
          cerr << "[adapters]\tadapter size is more then 32. Use slow adapters search" << "\n";
          do_adapter_optimized = false;
        }
      }

      // store information
      adapter_names.push_back(adapter_name);
      adapter_seqs.push_back(adapter_seq);
      adapter_hashes.push_back(hash_adapter(adapter_seq));

      if (adapter_size == 0) {
        adapter_size = adapter_seq.size();
        shortest_adapter_size = adapter_size;
      }
      else if (adapter_seq.size() != adapter_size) {
        cerr << "[adapters]\tadapters have different size. Use slow adapters search" << "\n";
        do_adapter_optimized = false;
        if(adapter_seq.size() < shortest_adapter_size){
          shortest_adapter_size = adapter_seq.size();
        }
      }

      line_by_space.clear();
    }
  }
  in.close();
}

void
FalcoConfig::read_contaminants_file() {
  if (contaminants_file == "") {
    if (!quiet)
      cerr << "[contaminants]\tusing default contaminant list (no file specified)\n";
    contaminants = FileConstants::contaminants;
    return;
  }
  ifstream in(contaminants_file);
  if (!in)
    throw runtime_error("contaminants file not found: " + contaminants_file);

  if (!quiet)
    cerr << "[contaminants]\tusing file " << contaminants_file << "\n";
  vector<string> line_by_space;

  // The contaminants file has a space separated name, and the last
  // instance is the biological sequence
  string line;
  contaminants.clear();
  while (getline(in, line)) {
    if (is_content_line(line)) {
      istringstream iss(line);
      string token;
      while (iss >> token)
        line_by_space.push_back(token);

      if (line_by_space.size() > 1) {
        string contaminant_name;
        for (size_t i = 0; i < line_by_space.size() - 1; ++i)
          contaminant_name += line_by_space[i] + " ";
        contaminants.push_back(make_pair(contaminant_name, line_by_space.back()));
      }
      line_by_space.clear();
    }
  }
}

const string FalcoConfig::html_template =
"<html>"
"<head>"
"    <meta charset=\"utf-8\">"
"    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1, shrink-to-fit=no\">"
""
"	<title>"
"     {{filename}} - report"
"	</title>"
"<link rel=\"stylesheet\" href=\"https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css\" integrity=\"sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T\" crossorigin=\"anonymous\">"
"<link href=\"https://stackpath.bootstrapcdn.com/font-awesome/4.7.0/css/font-awesome.min.css\" rel=\"stylesheet\" integrity=\"sha384-wvfXpqpZZVQGK6TAh5PVlGOfQNHSoD2xbE+QkPxCAFlNEevoEH3Sl0sibVcOQVnN\" crossorigin=\"anonymous\">"
"<style type=\"text/css\">"
" @media screen {"
"  div.summary {"
"    width: 18em;"
"    position:fixed;"
"    top: 4em;"
"    margin:1em 0 0 1em;"
"  }"
"  "
"  div.main {"
"    display:block;"
"    position:absolute;"
"    overflow:auto;"
"    height:auto;"
"    width:auto;"
"    top:4.5em;"
"    bottom:2.3em;"
"    left:18em;"
"    right:0;"
"    border-left: 1px solid #CCC;"
"    padding:0 0 0 1em;"
"    background-color: white;"
"    z-index:1;"
"  }"
"  "
"  div.header {"
"    background-color: #EEE;"
"    border:0;"
"    margin:0;"
"    padding: 0.2em;"
"    font-size: 200%;"
"    position:fixed;"
"    width:100%;"
"    top:0;"
"    left:0;"
"    z-index:2;"
"  }"
""
"  div.footer {"
"    background-color: #EEE;"
"    border:0;"
"    margin:0;"
"	padding:0.5em;"
"    height: 2.5em;"
"	overflow:hidden;"
"    font-size: 100%;"
"    position:fixed;"
"    bottom:0;"
"    width:100%;"
"    z-index:2;"
"  }"
"  "
"  img.indented {"
"    margin-left: 3em;"
"  }"
" }"
" "
" @media print {"
"	img {"
"		max-width:100% !important;"
"		page-break-inside: avoid;"
"	}"
"	h2, h3 {"
"		page-break-after: avoid;"
"	}"
"	div.header {"
"      background-color: #FFF;"
"    }"
"	"
" }"
" "
" body {    "
"  color: #000;   "
"  background-color: #FFF;"
"  border: 0;"
"  margin: 0;"
"  padding: 0;"
"  }"
"  "
"  div.header {"
"  border:0;"
"  margin:0;"
"  padding: 0.5em;"
"  font-size: 200%;"
"  width:100%;"
"  }    "
"  "
"  #header_title {"
"  display:inline-block;"
"  float:left;"
"  clear:left;"
"  }"
"  #header_filename {"
"  display:inline-block;"
"  float:right;"
"  clear:right;"
"  font-size: 50%;"
"  margin-right:2em;"
"  text-align: right;"
"  }"
""
"  div.header h3 {"
"  font-size: 50%;"
"  margin-bottom: 0;"
"  }"
"  "
"  div.summary ul {"
"  padding-left:0;"
"  list-style-type:none;"
"  }"
"  "
"  div.summary ul li img {"
"  margin-bottom:-0.5em;"
"  margin-top:0.5em;"
"  }"
"	  "
"  div.main {"
"  background-color: white;"
"  }"
"      "
"  div.module {"
"  padding-bottom:3em;"
"  padding-top:3em;"
"  border-bottom: 1px solid #990000"
"  }"
"	  "
"  div.footer {"
"  background-color: #EEE;"
"  border:0;"
"  margin:0;"
"  padding: 0.5em;"
"  font-size: 100%;"
"  width:100%;"
"  }"
""
"  h2 {"
"  color: #2a5e8c;"
"  padding-bottom: 0;"
"  margin-bottom: 0;"
"  clear:left;"
"  }"
""
"table {"
"  margin-left: 3em;"
"  text-align: center;"
"  }"
"  "
"th {"
"  text-align: center;"
"  background-color: #000080;"
"  color: #FFF;"
"  padding: 0.4em;"
"}"
"  "
"td {"
"  font-family: monospace;"
"  text-align: left;"
"  background-color: #EEE;"
"  color: #000;"
"  padding: 0.4em;"
"}"
""
"img {"
"  padding-top: 0;"
"  margin-top: 0;"
"  border-top: 0;"
"}"
""
"  "
"p {"
"  padding-top: 0;"
"  margin-top: 0;"
"}"
""
".pass {"
"  color : #009900;"
"}"
""
".warn {"
"  color : #999900;"
"}"
""
".fail {"
"  color : #990000;"
"}"
"</style>"
"<script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>"
"</head>"
"<body>"
"<div class=\"header\">"
"	<div id=\"header_title\">Report</div>"
"  <div id=\"header_filename\">{{date}}<br/> {{filename}}"
"	</div>"
"</div>"
"<div class=\"summary\"><h2>Summary</h2>"
"<ul>"
"  {{basicstatisticscs}}"
"  <li><a class=\"{{passbasicstatistics}}\" href=\"#basicstatistics\">"
"    {{basicstatisticsname}}"
"  </a></li>"
"  {{basicstatisticsce}}"
""
"  {{perbasesequencequalitycs}}"
"	<li><a class=\"{{passperbasesequencequality}}\" href=\"#perbasesequencequality\">"
"    {{perbasesequencequalityname}}</a></li>"
"  {{perbasesequencequalityce}}"
""
"  {{pertilesequencequalitycs}}"
"	<li><a class=\"{{passpertilesequencequality}}\" href=\"#pertilesequencequality\">{{pertilesequencequalityname}}</a></li>"
"  {{pertilesequencequalityce}}"
""
"  {{persequencequalityscorescs}}"
"	<li><a class=\"{{passpersequencequalityscores}}\" href=\"#persequencequalityscores\">{{persequencequalityscoresname}}</a></li>"
"  {{persequencequalityscoresce}}"
""
"  {{perbasesequencecontentcs}}"
"	<li><a class=\"{{passperbasesequencecontent}}\" href=\"#perbasesequencecontent\">{{perbasesequencecontentname}}</a></li>"
"  {{perbasesequencecontentce}}"
""
"  {{persequencegccontentcs}}"
"	<li><a class=\"{{passpersequencegccontent}}\" href=\"#persequencegccontent\">{{persequencegccontentname}}</a></li>"
"  {{persequencegccontentce}}"
""
"  {{perbasencontentcs}}"
"	<li><a class=\"{{passperbasencontent}}\" href=\"#perbasencontent\">{{perbasencontentname}}</a></li>"
"  {{perbasencontentce}}"
""
"  {{sequencelengthdistributioncs}}"
"	<li><a class=\"{{passsequencelengthdistribution}}\" href=\"#sequencelengthdistribution\">{{sequencelengthdistributionname}}</a></li>"
"  {{sequencelengthdistributionce}}"
""
"  {{sequenceduplicationlevelscs}}"
"	<li><a class=\"{{passsequenceduplicationlevels}}\" href=\"#sequenceduplicationlevels\">{{sequenceduplicationlevelsname}}</a></li>"
"  {{sequenceduplicationlevelsce}}"
""
"  {{overrepresentedsequencescs}}"
"	<li><a class=\"{{passoverrepresentedsequences}}\" href=\"#overrepresentedsequences\">{{overrepresentedsequencesname}}</a></li>"
"  {{overrepresentedsequencesce}}"
""
"  {{adaptercontentcs}}"
"	<li><a class=\"{{passadaptercontent}}\" href=\"#adaptercontent\">{{adaptercontentname}}</a></li>"
"  {{adaptercontentce}}"
""
"  {{kmercontentcs}}"
"	<li><a class=\"{{passkmercontent}}\" href=\"#kmercontent\">{{kmercontentname}}</a></li>"
"  {{kmercontentce}}"
""
""
"</ul>"
"</div>"
"<div class=\"main\">"
"<div class=\"module\">"
"  <h2 class=\"{{passbasicstatistics}}\" id=\"basicstatistics\">"
"    {{basicstatisticsname}}: {{passbasicstatistics}}"
"  </h2>"
"  {{basicstatisticsdata}}"
"</div>"
""
"{{perbasesequencequalitycs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passperbasesequencequality}}\" id=\"perbasesequencequality\">"
"    {{perbasesequencequalityname}}: {{passperbasesequencequality}}</h2>"
" 	<div id=\"seqbasequalityboxplot\"></div>"
"</div>"
"{{perbasesequencequalityce}}"
""
""
"{{pertilesequencequalitycs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passpertilesequencequality}}\" id=\"pertilesequencequality\">"
"    {{pertilesequencequalityname}} : {{passpertilesequencequality}}"
"  </h2>"
" 	<div id=\"tilequalityheatmap\"></div>"
"</div>"
"{{pertilesequencequalityce}}"
""
"{{persequencequalityscorescs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passpersequencequalityscores}}\" id=\"persequencequalityscores\">"
"    {{persequencequalityscoresname}} : {{passpersequencequalityscores}}"
"  </h2>"
" 	<div id=\"seqqualitylineplot\"></div>"
"</div>"
"{{persequencequalityscoresce}}"
""
"{{perbasesequencecontentcs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passperbasesequencecontent}}\" id=\"perbasesequencecontent\">"
"    {{perbasesequencecontentname}} : {{passperbasesequencecontent}}"
"  </h2>"
" 	<div id=\"basesequencecontentlineplot\"></div>"
"</div>"
"{{perbasesequencecontentce}}"
""
"{{persequencegccontentcs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passpersequencegccontent}}\" id=\"persequencegccontent\">"
"    {{persequencegccontentname}}: {{passpersequencegccontent}}"
"  </h2>"
" 	<div id=\"sequencegccontentlineplot\"></div>"
"</div>"
"{{persequencegccontentce}}"
""
""
"{{perbasencontentcs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passperbasencontent}}\" id=\"perbasencontent\">"
"    {{perbasencontentname}} : {{passperbasencontent}}"
"  </h2>"
" 	<div id=\"basencontentlineplot\"></div>"
"</div>"
"{{perbasencontentce}}"
""
"{{sequencelengthdistributioncs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passsequencelengthdistribution}}\" id=\"sequencelengthdistribution\">"
"    {{sequencelengthdistributionname}} : {{passsequencelengthdistribution}}"
"  </h2>"
" 	<div id=\"sequencelengthdistributionlineplot\"></div>"
"</div>"
"{{sequencelengthdistributionce}}"
""
"{{sequenceduplicationlevelscs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passsequenceduplicationlevels}}\" id=\"sequenceduplicationlevels\">"
"    {{sequenceduplicationlevelsname}} : {{passsequenceduplicationlevels}}"
"  </h2>"
" 	<div id=\"seqduplevelslineplot\"></div>"
"</div>"
"{{sequenceduplicationlevelsce}}"
""
"{{overrepresentedsequencescs}}"
"<div class=\"module\">"
"	<h2 class=\"{{passoverrepresentedsequences}}\" id=\"overrepresentedsequences\">"
"    {{overrepresentedsequencesname}} : {{passoverrepresentedsequences}}</h2>"
"  {{overrepresentedsequencesdata}}"
"</div>"
"{{overrepresentedsequencesce}}"
""
"{{adaptercontentcs}}"
"<div class=\"module\">"
"  <h2 class=\"{{passadaptercontent}}\" id=\"adaptercontent\">"
"    {{adaptercontentname}} : {{passadaptercontent}}"
"  </h2>"
" 	<div id=\"adapterlineplot\"></div>"
"</div>"
"{{adaptercontentce}}"
""
"{{kmercontentcs}}"
"<div class=\"module\">"
"  <h2 class=\"{{passkmercontent}}\" id=\"kmercontent\">"
"    {{kmercontentname}} : {{passkmercontent}}"
"  </h2>"
" 	<div id=\"kmerlineplot\"></div>"
"</div>"
"{{kmercontentce}}"
""
""
"</div>"
"<div class=\"footer\">Falco " + FalcoConfig::FalcoVersion +
"</div>"
"</body>"
"<script src=\"https://code.jquery.com/jquery-3.3.1.slim.min.js\" "
"integrity=\"sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo\" "
"crossorigin=\"anonymous\"></script>"
""
"<script src=\"https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.7/umd/popper.min.js\" "
"integrity=\"sha384-UO2eT0CpHqdSJQ6hJty5KVphtPhzWj9WO1clHTMGa3JDZwrnQq4sF86dIHNDz0W1\""
"crossorigin=\"anonymous\"></script>"
""
"<script src=\"https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/js/bootstrap.min.js\""
"integrity=\"sha384-JjSmVgyd0p3pXB1rRibZUAYoIIy6OrQ6VrjIEaFf/nJGzIxFDsf4x0xIM+B07jRM\""
"crossorigin=\"anonymous\"></script>"
""
"<script>"
"  "
"if (document.getElementById('seqbasequalityboxplot') !== null) {"
"  Plotly.newPlot('seqbasequalityboxplot', ["
"   {{perbasesequencequalitydata}}"
"  ], {"
"    margin: { t: 0 }, showlegend: false,"
"    xaxis : {title : 'Base position'},"
"    yaxis : {title : 'Phread quality'},"
"  });"
"}"
"if (document.getElementById('tilequalityheatmap') !== null) {"
"  Plotly.newPlot('tilequalityheatmap', ["
"   {{pertilesequencequalitydata}}"
"  ], {"
"    margin: { t: 0 }, "
"    showlegend: false,"
"    xaxis : {title : 'Base position'},"
"    yaxis : {title : 'tile', type: 'category'} "
"  });"
"}"
""
"if (document.getElementById('seqqualitylineplot') !== null) {"
"  Plotly.newPlot('seqqualitylineplot', ["
"   {{persequencequalityscoresdata}}"
"  ], {"
"    margin: { t: 0 },"
"    showlegend: true,"
"    xaxis : {title : 'Phread quality'},"
"    yaxis : {title : 'Density'}"
"  } );"
"}"
"if (document.getElementById('basesequencecontentlineplot') !== null) {"
"  Plotly.newPlot('basesequencecontentlineplot', ["
"   {{perbasesequencecontentdata}}"
"  ], {"
"    margin: { t: 0 },"
"    showlegend: true,"
"    xaxis : {title : 'Base position'},"
"    yaxis : {title : '% sequence content'}"
"  } );"
"}"
"if (document.getElementById('sequencegccontentlineplot') !== null) {"
"  Plotly.newPlot('sequencegccontentlineplot', ["
"   {{persequencegccontentdata}}"
"  ], {"
"    margin: { t: 0 }, "
"    showlegend: true,"
"    xaxis : {title : '% GC'},"
"    yaxis : {title : 'Density'}"
"  } );"
"}"
"if (document.getElementById('basencontentlineplot') !== null) {"
"  Plotly.newPlot('basencontentlineplot', ["
"   {{perbasencontentdata}}"
"  ], {"
"    margin: { t: 0 }, "
"    showlegend: true,"
"    xaxis : {title : 'Base position'},"
"    yaxis : {title : '% N'}"
"  } );"
"}"
""
"if (document.getElementById('sequencelengthdistributionlineplot') !== null) {"
"  Plotly.newPlot('sequencelengthdistributionlineplot', ["
"   {{sequencelengthdistributiondata}}"
"  ], {"
"    margin: { t: 0 }, "
"    showlegend: true,"
"    xaxis : {title : 'Sequence length'},"
"    yaxis : {title : 'Number of sequences'}"
"  } );"
"}"
"if (document.getElementById('seqduplevelslineplot') !== null) {"
"  Plotly.newPlot('seqduplevelslineplot', ["
"   {{sequenceduplicationlevelsdata}}"
"  ], {"
"    margin: { t: 0 }, "
"    showlegend: true,"
"    xaxis : {title : 'Duplication rate',"
"             tickvals : [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16],"
"             ticktext : ['1','2','3','4','5','6','7','8','9','10+','50+','100+','500+','1k+','5k+','10k+']},"
"    yaxis : {title : '% of sequences'}"
"  } );"
"}"
""
"if (document.getElementById('adapterlineplot') !== null) {"
"  Plotly.newPlot('adapterlineplot', ["
"   {{adaptercontentdata}}"
"  ], {"
"    margin: { t: 0 }, "
"    showlegend: true,"
"    xaxis : {title : 'Base position'},"
"    yaxis : {title : '% sequences with adapter before position'}"
"  } );"
"}"
"if (document.getElementById('kmerlineplot') !== null) {"
"  Plotly.newPlot('kmerlineplot', ["
"   {{kmercontentdata}}"
"  ], {"
"    margin: { t: 0 }, "
"    showlegend: true,"
"    xaxis : {title : 'Base position'},"
"    yaxis : {title : 'log2(obs/ exp max)'}"
"  } );"
"}"

"</script>"
"</html>";
