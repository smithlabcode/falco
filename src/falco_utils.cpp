// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "falco_utils.hpp"

#include <cmath>
#include <cstdint>
#include <ctime>
#include <format>
#include <iomanip>
#include <sstream>
#include <string>

[[nodiscard]] auto
size_to_units(const std::int64_t s, const std::string &suffix) -> std::string {
  const auto as_frac_2 = [](const auto a, const auto b) {
    return std::floor(10 * as_frac(a, b)) / 10;
  };
  if (s >= gigabytes)
    return std::format("{}G{}", as_frac_2(s, gigabytes), suffix);
  if (s >= megabytes)
    return std::format("{}M{}", as_frac_2(s, megabytes), suffix);
  if (s >= kilobytes)
    return std::format("{}K{}", as_frac_2(s, kilobytes), suffix);
  return std::format("{}", s);
};

[[nodiscard]] auto
get_program_start_time() -> std::chrono::time_point<std::chrono::system_clock> {
  static const auto start_time = std::chrono::system_clock::now();
  return start_time;
}

[[nodiscard]] auto
format_program_start_date_and_time() -> std::string {
  const auto t = get_program_start_time();
  const auto t_c = std::chrono::system_clock::to_time_t(t);
  std::ostringstream oss;
  // NOLINTNEXTLINE(concurrency-mt-unsafe)
  oss << std::put_time(std::localtime(&t_c), "%F %T %Z");
  return oss.str();
}
