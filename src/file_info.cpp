// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "file_info.hpp"
#include "falco_utils.hpp"

#include <format>
#include <string>

[[nodiscard]] auto
file_info::to_line() const -> std::string {
  return std::format("{}\t{}\t{}", name,
                     description + (!has_quals ? std::string(" (quals missing)")
                                               : std::string{}),
                     size_to_units(size, std::string{}));
}
