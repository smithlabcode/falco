// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "contaminant_set.hpp"
#include "falco_utils.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iterator>
#include <ranges>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>  // for std::get
#include <utility>
#include <vector>

// clang-format off
static constexpr auto default_contaminants = {std::pair
  {R"(Illumina Single End Adapter 1)", R"(GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(Illumina Single End Adapter 2)", R"(CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT)"},
  {R"(Illumina Single End PCR Primer 1)", R"(AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Single End PCR Primer 2)", R"(CAAGCAGAAGACGGCATACGAGCTCTTCCGATCT)"},
  {R"(Illumina Single End Sequencing Primer)", R"(ACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Paired End Adapter 1)", R"(ACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Paired End Adapter 2)", R"(GATCGGAAGAGCGGTTCAGCAGGAATGCCGAG)"},
  {R"(Illumina Paried End PCR Primer 1)", R"(AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Paired End PCR Primer 2)", R"(CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT)"},
  {R"(Illumina Paried End Sequencing Primer 1)", R"(ACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Paired End Sequencing Primer 2)", R"(CGGTCTCGGCATTCCTGCTGAACCGCTCTTCCGATCT)"},
  {R"(Illumina DpnII expression Adapter 1)", R"(ACAGGTTCAGAGTTCTACAGTCCGAC)"},
  {R"(Illumina DpnII expression Adapter 2)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina DpnII expression PCR Primer 1)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina DpnII expression PCR Primer 2)", R"(AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA)"},
  {R"(Illumina DpnII expression Sequencing Primer)", R"(CGACAGGTTCAGAGTTCTACAGTCCGACGATC)"},
  {R"(Illumina NlaIII expression Adapter 1)", R"(ACAGGTTCAGAGTTCTACAGTCCGACATG)"},
  {R"(Illumina NlaIII expression Adapter 2)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina NlaIII expression PCR Primer 1)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina NlaIII expression PCR Primer 2)", R"(AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA)"},
  {R"(Illumina NlaIII expression Sequencing Primer)", R"(CCGACAGGTTCAGAGTTCTACAGTCCGACATG)"},
  {R"(Illumina Small RNA Adapter 1)", R"(GTTCAGAGTTCTACAGTCCGACGATC)"},
  {R"(Illumina Small RNA Adapter 2)", R"(TGGAATTCTCGGGTGCCAAGG)"},
  {R"(Illumina Small RNA RT Primer)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina Small RNA PCR Primer 2)", R"(AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA)"},
  {R"(Illumina Small RNA Sequencing Primer)", R"(CGACAGGTTCAGAGTTCTACAGTCCGACGATC)"},
  {R"(Illumina Multiplexing Adapter 1)", R"(GATCGGAAGAGCACACGTCT)"},
  {R"(Illumina Multiplexing Adapter 2)", R"(ACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Multiplexing PCR Primer 1.01)", R"(AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Multiplexing PCR Primer 2.01)", R"(GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT)"},
  {R"(Illumina Multiplexing Read1 Sequencing Primer)", R"(ACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(Illumina Multiplexing Index Sequencing Primer)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCAC)"},
  {R"(Illumina Multiplexing Read2 Sequencing Primer)", R"(GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT)"},
  {R"(Illumina PCR Primer Index 1)", R"(CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 2)", R"(CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 3)", R"(CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 4)", R"(CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 5)", R"(CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 6)", R"(CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 7)", R"(CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 8)", R"(CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 9)", R"(CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 10)", R"(CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 11)", R"(CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTC)"},
  {R"(Illumina PCR Primer Index 12)", R"(CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTC)"},
  {R"(Illumina DpnII Gex Adapter 1)", R"(GATCGTCGGACTGTAGAACTCTGAAC)"},
  {R"(Illumina DpnII Gex Adapter 1.01)", R"(ACAGGTTCAGAGTTCTACAGTCCGAC)"},
  {R"(Illumina DpnII Gex Adapter 2)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina DpnII Gex Adapter 2.01)", R"(TCGTATGCCGTCTTCTGCTTG)"},
  {R"(Illumina DpnII Gex PCR Primer 1)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina DpnII Gex PCR Primer 2)", R"(AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA)"},
  {R"(Illumina DpnII Gex Sequencing Primer)", R"(CGACAGGTTCAGAGTTCTACAGTCCGACGATC)"},
  {R"(Illumina NlaIII Gex Adapter 1.01)", R"(TCGGACTGTAGAACTCTGAAC)"},
  {R"(Illumina NlaIII Gex Adapter 1.02)", R"(ACAGGTTCAGAGTTCTACAGTCCGACATG)"},
  {R"(Illumina NlaIII Gex Adapter 2.01)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina NlaIII Gex Adapter 2.02)", R"(TCGTATGCCGTCTTCTGCTTG)"},
  {R"(Illumina NlaIII Gex PCR Primer 1)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(Illumina NlaIII Gex PCR Primer 2)", R"(AATGATACGGCGACCACCGACAGGTTCAGAGTTCTACAGTCCGA)"},
  {R"(Illumina NlaIII Gex Sequencing Primer)", R"(CCGACAGGTTCAGAGTTCTACAGTCCGACATG)"},
  {R"(Illumina 5p RNA Adapter)", R"(GTTCAGAGTTCTACAGTCCGACGATC)"},
  {R"(Illumina RNA Adapter1)", R"(TGGAATTCTCGGGTGCCAAGG)"},
  {R"(Illumina Small RNA 3p Adapter 1)", R"(ATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(Illumina Small RNA PCR Primer 1)", R"(CAAGCAGAAGACGGCATACGA)"},
  {R"(TruSeq Universal Adapter)", R"(AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT)"},
  {R"(TruSeq Adapter, Index 1)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 2)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATGTATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 3)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACTTAGGCATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 4)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 5)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACACAGTGATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 6)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCAATATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 7)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACCAGATCATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 8)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTTGAATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 9)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACGATCAGATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 10)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACTAGCTTATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 11)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGCTACATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 12)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACCTTGTAATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 13)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 14)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTTCCGTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 15)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACATGTCAGTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 16)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCGTCCCTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 18)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 19)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGAAACTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 20)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTGGCCTTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 21)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACGTTTCGGTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 22)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACCGTACGTTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 23)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACCCACTCTTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 25)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACACTGATATCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(TruSeq Adapter, Index 27)", R"(GATCGGAAGAGCACACGTCTGAACTCCAGTCACATTCCTTTCTCGTATGCCGTCTTCTGCTTG)"},
  {R"(Illumina RNA RT Primer)", R"(GCCTTGGCACCCGAGAATTCCA)"},
  {R"(Illumina RNA PCR Primer)", R"(AATGATACGGCGACCACCGAGATCTACACGTTCAGAGTTCTACAGTCCGA)"},
  {R"(RNA PCR Primer, Index 1)", R"(CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 2)", R"(CAAGCAGAAGACGGCATACGAGATACATCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 3)", R"(CAAGCAGAAGACGGCATACGAGATGCCTAAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 4)", R"(CAAGCAGAAGACGGCATACGAGATTGGTCAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 5)", R"(CAAGCAGAAGACGGCATACGAGATCACTGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 6)", R"(CAAGCAGAAGACGGCATACGAGATATTGGCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 7)", R"(CAAGCAGAAGACGGCATACGAGATGATCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 8)", R"(CAAGCAGAAGACGGCATACGAGATTCAAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 9)", R"(CAAGCAGAAGACGGCATACGAGATCTGATCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 10)", R"(CAAGCAGAAGACGGCATACGAGATAAGCTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 11)", R"(CAAGCAGAAGACGGCATACGAGATGTAGCCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 12)", R"(CAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 13)", R"(CAAGCAGAAGACGGCATACGAGATTTGACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 14)", R"(CAAGCAGAAGACGGCATACGAGATGGAACTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 15)", R"(CAAGCAGAAGACGGCATACGAGATTGACATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 16)", R"(CAAGCAGAAGACGGCATACGAGATGGACGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 17)", R"(CAAGCAGAAGACGGCATACGAGATCTCTACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 18)", R"(CAAGCAGAAGACGGCATACGAGATGCGGACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 19)", R"(CAAGCAGAAGACGGCATACGAGATTTTCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 20)", R"(CAAGCAGAAGACGGCATACGAGATGGCCACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 21)", R"(CAAGCAGAAGACGGCATACGAGATCGAAACGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 22)", R"(CAAGCAGAAGACGGCATACGAGATCGTACGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 23)", R"(CAAGCAGAAGACGGCATACGAGATCCACTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 24)", R"(CAAGCAGAAGACGGCATACGAGATGCTACCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 25)", R"(CAAGCAGAAGACGGCATACGAGATATCAGTGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 26)", R"(CAAGCAGAAGACGGCATACGAGATGCTCATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 27)", R"(CAAGCAGAAGACGGCATACGAGATAGGAATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 28)", R"(CAAGCAGAAGACGGCATACGAGATCTTTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 29)", R"(CAAGCAGAAGACGGCATACGAGATTAGTTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 30)", R"(CAAGCAGAAGACGGCATACGAGATCCGGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 31)", R"(CAAGCAGAAGACGGCATACGAGATATCGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 32)", R"(CAAGCAGAAGACGGCATACGAGATTGAGTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 33)", R"(CAAGCAGAAGACGGCATACGAGATCGCCTGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 34)", R"(CAAGCAGAAGACGGCATACGAGATGCCATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 35)", R"(CAAGCAGAAGACGGCATACGAGATAAAATGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 36)", R"(CAAGCAGAAGACGGCATACGAGATTGTTGGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 37)", R"(CAAGCAGAAGACGGCATACGAGATATTCCGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 38)", R"(CAAGCAGAAGACGGCATACGAGATAGCTAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 39)", R"(CAAGCAGAAGACGGCATACGAGATGTATAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 40)", R"(CAAGCAGAAGACGGCATACGAGATTCTGAGGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 41)", R"(CAAGCAGAAGACGGCATACGAGATGTCGTCGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 42)", R"(CAAGCAGAAGACGGCATACGAGATCGATTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 43)", R"(CAAGCAGAAGACGGCATACGAGATGCTGTAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 44)", R"(CAAGCAGAAGACGGCATACGAGATATTATAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 45)", R"(CAAGCAGAAGACGGCATACGAGATGAATGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 46)", R"(CAAGCAGAAGACGGCATACGAGATTCGGGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 47)", R"(CAAGCAGAAGACGGCATACGAGATCTTCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(RNA PCR Primer, Index 48)", R"(CAAGCAGAAGACGGCATACGAGATTGCCGAGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA)"},
  {R"(ABI Dynabead EcoP Oligo)", R"(CTGATCTAGAGGTACCGGATCCCAGCAGT)"},
  {R"(ABI Solid3 Adapter A)", R"(CTGCCCCGGGTTCCTCATTCTCTCAGCAGCATG)"},
  {R"(ABI Solid3 Adapter B)", R"(CCACTACGCCTCCGCTTTCCTCTCTATGGGCAGTCGGTGAT)"},
  {R"(ABI Solid3 5' AMP Primer)", R"(CCACTACGCCTCCGCTTTCCTCTCTATG)"},
  {R"(ABI Solid3 3' AMP Primer)", R"(CTGCCCCGGGTTCCTCATTCT)"},
  {R"(ABI Solid3 EF1 alpha Sense Primer)", R"(CATGTGTGTTGAGAGCTTC)"},
  {R"(ABI Solid3 EF1 alpha Antisense Primer)", R"(GAAAACCAAAGTGGTCCAC)"},
  {R"(ABI Solid3 GAPDH Forward Primer)", R"(TTAGCACCCCTGGCCAAGG)"},
  {R"(ABI Solid3 GAPDH Reverse Primer)", R"(CTTACTCCTTGGAGGCCATG)"},
  {R"(Clontech Universal Primer Mix Short)", R"(CTAATACGACTCACTATAGGGC)"},
  {R"(Clontech Universal Primer Mix Long)", R"(CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT)"},
  {R"(Clontech SMARTer II A Oligonucleotide)", R"(AAGCAGTGGTATCAACGCAGAGTAC)"},
  {R"(Clontech SMART CDS Primer II A)", R"(AAGCAGTGGTATCAACGCAGAGTACT)"},
  {R"(Clontech_Universal_Primer_Mix_Short)", R"(CTAATACGACTCACTATAGGGC)"},
  {R"(Clontech_Universal_Primer_Mix_Long)", R"(CTAATACGACTCACTATAGGGCAAGCAGTGGTATCAACGCAGAGT)"},
  {R"(Clontech_SMARTer_II_A_Oligonucleotide)", R"(AAGCAGTGGTATCAACGCAGAGTAC)"},
  {R"(Clontech_SMART_CDS_Primer_II_A)", R"(AAGCAGTGGTATCAACGCAGAGTACT)"},
  {R"(Clontech_SMART_CDS_Primer_II_A)", R"(ACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGTCGC)"},
  {R"(Clontech_SMART_CDS_Primer_II_A)", R"(GAGTACGTACTCTGCGTTGATACCACTGCTTCCGCGGACAGGCGTGTAGATCTCGGTGGT)"},
};
// clang-format on

[[nodiscard]] static auto
load_contaminants(const std::string &filename)
  -> std::vector<std::pair<std::string, std::string>> {
  // ADS: (todo) handle carriage returns and other control chars
  std::ifstream in(filename);
  if (!in)
    throw std::runtime_error("failed to open contaminants file: " + filename);
  std::vector<std::pair<std::string, std::string>> contaminants;
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
    for (auto itr = std::cbegin(line); itr != std::cend(line);
         itr = std::next(itr))
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
  return contaminants;
}

[[nodiscard]] auto
contaminant_set::get_name(const std::int64_t contam_idx)
  -> const std::string & {
  using std::string_literals::operator""s;
  static constexpr auto no_hit_label = "No Hit"s;
  const auto &contaminants = instance().contaminants;
  if (contam_idx < 0 || contam_idx >= std::ssize(contaminants))
    return no_hit_label;
  return contaminants[contam_idx].first;
}

/// get the longest substring of left that is a prefix of right
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
contaminant_set::match(const std::string &query) -> std::int64_t {
  const auto &contaminants = instance().contaminants;
  auto best_idx = 0L;
  auto best_match = 0L;
  auto best_match_len = 0L;
  for (const auto &[idx, seq] :
       falco::views::enumerate(std::views::elements<1>(contaminants))) {
    const auto n_match =
      std::max(get_overlap(query, seq), get_overlap(seq, query));
    if (n_match > best_match) {
      best_idx = idx;
      best_match = n_match;
      best_match_len = std::ssize(seq);
    }
  }
  const auto match_cutoff = std::min(best_match_len, std::ssize(query)) / 2;
  // If any sequence is a match, return the best one
  return best_match < match_cutoff ? -1 : best_idx;
}

contaminant_set::contaminant_set(const std::string &filename) {
  if (filename.empty())
    std::ranges::copy(default_contaminants, std::back_inserter(contaminants));
  else
    contaminants = load_contaminants(filename);
}
