// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_ORIGINAL_DUPLICATES_HPP_
#define SRC_ORIGINAL_DUPLICATES_HPP_

#include "duplication_results.hpp"
#include "falco_word.hpp"
#include "file_info.hpp"

#include <string>
#include <vector>

[[nodiscard]] auto
initialize_original_duplicates(
  const std::vector<std::string> &infiles, const std::vector<file_info> &infos,
  const std::uint32_t n_threads) -> std::vector<dups_init_t>;

#endif  // SRC_ORIGINAL_DUPLICATES_HPP_
