// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_CONTAMINANT_SET_HPP_
#define SRC_CONTAMINANT_SET_HPP_

#include <cstdint>
#include <iterator>
#include <ranges>  // for std::pair
#include <string>
#include <utility>
#include <vector>

struct contaminant_set {
  static auto
  instance(const std::string &filename = std::string{})
    -> const contaminant_set & {
    static const contaminant_set s(filename);
    return s;
  }

  [[nodiscard]] static auto
  n_contaminants() -> std::uint64_t {
    return std::size(instance().contaminants);
  }

  [[nodiscard]] static auto
  get_name(std::int64_t idx) -> const std::string &;

  [[nodiscard]] static auto
  match(const std::string &query) -> std::int64_t;

  // clang-format off
  contaminant_set(const contaminant_set &) = delete;
  contaminant_set(contaminant_set &&) = delete;
  auto operator=(const contaminant_set &) = delete;
  auto operator=(contaminant_set &&) = delete;
  ~contaminant_set() = default;
  // clang-format on

private:
  std::vector<std::pair<std::string, std::string>> contaminants;
  explicit contaminant_set(const std::string &filename);
};  // contaminant_set

#endif  // SRC_CONTAMINANT_SET_HPP_
