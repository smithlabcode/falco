/* MIT License
 *
 * Copyright (c) 2026 Andrew D Smith
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef SRC_ADAPTER_SET_HPP_
#define SRC_ADAPTER_SET_HPP_

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
  instance(const std::string &filename = std::string{}) -> const adapter_set & {
    static const adapter_set s(filename);
    return s;
  }

  [[nodiscard]] auto
  validate() const -> std::tuple<bool, std::string>;

  adapter_set(const adapter_set &) = delete;
  auto
  operator=(const adapter_set &) -> adapter_set & = delete;

private:
  explicit adapter_set(const std::string &filename);
  adapter_set() = default;
  ~adapter_set() = default;
};  // adapter_set

#endif  // SRC_ADAPTER_SET_HPP_
