// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "run_mode.hpp"

#include <algorithm>
#include <format>
#include <iterator>
#include <ranges>
#include <string>
#include <unordered_map>
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
//
// ADS (2026-08-15 edit): additional modes have been added, but I'm not going to
// adjust the above right now because of time. So here are the extras, none of
// which can be set in the config file, only cli:
//
// do_groups (had been part of run_mode)
// do_bisulfite (changes grading and nothing else)
// do_original_dups (use the original duplication mode from FastQC and Falco v1)
// do_preseq (make another output file for input to preseq)
//
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
run_mode::assign(const std::unordered_map<std::string, bool> &modes) -> void {
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
  // settings below are not in config file
  if (do_groups_ == 0) do_groups_ = do_groups_default;
  if (do_bisulfite_ == 0) do_bisulfite_ = do_bisulfite_default;
  if (do_preseq_ == 0) do_preseq_ = do_preseq_default;
  if (do_original_dups_ == 0) do_original_dups_ = do_original_dups_default;
  //
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
  // clang-format on
}

[[nodiscard]] auto
run_mode::string_verbose() const -> std::string {
  static constexpr auto spacer = 2;
  // ADS: need a mechanism to keep this map it in sync with var names and labels
  static const auto descriptions = std::unordered_map<std::string, std::string>{
    {"adapter", "adapter content"},
    {"duplication", "sequence duplication"},
    {"gc_sequence", "GC content"},
    {"kmer", "k-mer content"},
    {"n_content", "N content"},
    {"overrepresented", "overrep sequences"},
    {"quality_base", "per-base quality"},
    {"quality_sequence", "per-sequence quality"},
    {"sequence", "sequence composition"},
    {"sequence_length", "sequence length"},
    {"tile", "per-tile quality"},
  };
  static const auto fmt_mode = [&](const std::string &label,
                                   const auto the_mode) {
    const auto itr = descriptions.find(label);
    if (itr == std::cend(descriptions))
      throw std::runtime_error("failed to find mode: " + label);
    return std::pair{itr->second, std::string{the_mode ? "ON" : "OFF"}};
  };
  std::vector<std::pair<std::string, std::string>> lines;
  lines.push_back(fmt_mode("adapter", do_adap()));
  lines.push_back(fmt_mode("duplication", do_dups()));
  lines.push_back(fmt_mode("gc_sequence", do_gc_content()));
  lines.push_back(fmt_mode("kmer", do_kmers()));
  lines.push_back(fmt_mode("sequence_length", do_length()));
  lines.push_back(fmt_mode("n_content", do_n_content()));
  lines.push_back(fmt_mode("overrepresented", do_overrep()));
  lines.push_back(fmt_mode("quality_base", do_qual_base()));
  lines.push_back(fmt_mode("quality_sequence", do_qual_seq()));
  lines.push_back(fmt_mode("sequence", do_sequence()));
  lines.push_back(fmt_mode("tile", do_tiles()));
  const auto w = spacer + std::ranges::max(std::views::transform(
                            std::views::values(descriptions),
                            [](const auto &x) { return std::size(x); }));
  std::string r;
  for (const auto &[descr, val] : lines)
    r += std::format("{:{}} {}\n", std::format("{}:", descr), w, val);
  return r;
}
