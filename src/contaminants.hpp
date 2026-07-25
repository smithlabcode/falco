// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_CONTAMINANTS_HPP_
#define SRC_CONTAMINANTS_HPP_

#include <cstdint>
#include <string>
#include <utility>  // IWYU pragma: keep
#include <vector>

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
extern std::vector<std::pair<std::string, std::string>> contaminants;

[[nodiscard]] auto
get_contam_name(std::int64_t contam_idx) -> const std::string &;

auto
load_contaminants(const std::string &filename) -> void;

[[nodiscard]] auto
match_contaminant(const std::string &query) -> std::int64_t;

#endif  // SRC_CONTAMINANTS_HPP_
