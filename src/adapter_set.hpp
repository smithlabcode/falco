// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_ADAPTER_SET_HPP_
#define SRC_ADAPTER_SET_HPP_

#include "run_mode.hpp"

#include <array>
#include <iterator>
#include <string>
#include <tuple>
#include <variant>
#include <vector>

struct adapter_set {
  static constexpr auto min_adapter_size = 6;
  static constexpr auto max_adapter_size = 16;

  static constexpr auto default_n_adapters = 6;
  static constexpr auto default_adapter_size = 12;
  static constexpr auto default_adapter_names = std::array{
    // clang-format off
    "Illumina Universal Adapter",
    "Illumina Small RNA 3' Adapter",
    "Illumina Small RNA 5' Adapter",
    "Nextera Transposase Sequence",
    "PolyA",
    "PolyG",
    // clang-format on
  };
  static constexpr auto default_adapters = std::array{
    // clang-format off
    "AGATCGGAAGAG",  // Illumina Universal Adapter
    "TGGAATTCTCGG",  // Illumina Small RNA 3' Adapter
    "GATCGTCGGACT",  // Illumina Small RNA 5' Adapter
    "CTGTCTCTTATA",  // Nextera Transposase Sequence
    "AAAAAAAAAAAA",  // PolyA
    "GGGGGGGGGGGG",  // PolyG
    // clang-format on
  };

  std::vector<std::string> adapter_names;
  std::vector<std::string> adapters;

  [[nodiscard]] static auto
  n_adapters() {
    return std::size(instance().adapters);
  }

  [[nodiscard]] static auto
  adapter_size() {
    return std::size(instance().adapters.front());
  }

  static auto
  instance(const run_mode &mode = {},
           const std::string &filename = std::string{}) -> const adapter_set & {
    static const adapter_set s(mode, filename);
    return s;
  }

  [[nodiscard]] auto
  validate() const -> std::tuple<bool, std::string>;

  adapter_set(const adapter_set &) = delete;
  adapter_set(adapter_set &&) = delete;
  auto
  operator=(const adapter_set &) = delete;
  auto
  operator=(adapter_set &&) = delete;
  ~adapter_set() = default;

private:
  explicit adapter_set(const run_mode &mode, const std::string &filename);
};  // adapter_set

#endif  // SRC_ADAPTER_SET_HPP_
