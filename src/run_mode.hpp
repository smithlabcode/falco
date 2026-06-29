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

#ifndef SRC_RUN_MODE_HPP_
#define SRC_RUN_MODE_HPP_

#include "boost/boost_unordered.hpp"

#include <string>
#include <vector>

class run_mode {
public:
  auto
  set_unassigned() -> void;

  [[nodiscard]] auto
  string() const -> std::string {
    return {};
  }

  // clang-format off
  // ADS: do_groups is not set in config file
  [[nodiscard]] auto do_groups() const -> bool { return do_groups_ == 1; }
  //
  [[nodiscard]] auto do_adap() const -> bool { return do_adap_ == 1; }
  [[nodiscard]] auto do_dups() const -> bool { return do_dups_ == 1; }
  [[nodiscard]] auto do_gc_content() const -> bool { return do_gc_content_ == 1; }
  [[nodiscard]] auto do_kmers() const -> bool { return do_kmers_ == 1; }
  [[nodiscard]] auto do_length() const -> bool { return do_length_ == 1; }
  [[nodiscard]] auto do_n_content() const -> bool { return do_n_content_ == 1; }
  [[nodiscard]] auto do_overrep() const -> bool { return do_overrep_ == 1; }
  [[nodiscard]] auto do_qual_base() const -> bool { return do_qual_base_ == 1; }
  [[nodiscard]] auto do_qual_seq() const -> bool { return do_qual_seq_ == 1; }
  [[nodiscard]] auto do_sequence() const -> bool { return do_sequence_ == 1; }
  [[nodiscard]] auto do_tiles() const -> bool { return do_tiles_ == 1; }
  // ADS: below, run duplication analysis if either dups or overrep is requested
  [[nodiscard]] auto do_dup_analysis() const -> bool {
    return do_overrep() || do_dups();
  }
  // clang-format on

  // clang-format off
  auto set_do_groups(const int x) { if (x) do_groups_ = x; }
  //
  auto set_do_adap(const int x) { if (x) do_adap_ = x; }
  auto set_do_dups(const int x) { if (x) do_dups_ = x; }
  auto set_do_gc_content(const int x) { if (x) do_gc_content_ = x; }
  auto set_do_kmers(const int x) { if (x) do_kmers_ = x; }
  auto set_do_length(const int x) { if (x) do_length_ = x; }
  auto set_do_n_content(const int x) { if (x) do_n_content_ = x; }
  auto set_do_overrep(const int x) { if (x) do_overrep_ = x; }
  auto set_do_qual_base(const int x) { if (x) do_qual_base_ = x; }
  auto set_do_qual_seq(const int x) { if (x) do_qual_seq_ = x; }
  auto set_do_sequence(const int x) { if (x) do_sequence_ = x; }
  auto set_do_tiles(const int x) { if (x) do_tiles_ = x; }

  [[nodiscard]] static auto
  get_labels() -> const std::vector<std::string> & { return labels; }
  // clang-format on

  auto
  assign(const boost::unordered_flat_map<std::string, bool> &modes) -> void;

private:
  // ADS: 1 is yes; -1 is no; 0 is not assigned
  static constexpr auto do_groups_default = -1;  // OFF

  static constexpr auto do_adap_default = 1;  // affects processing
  static constexpr auto do_dups_default = 1;  // affects processing
  static constexpr auto do_gc_content_default = 1;
  static constexpr auto do_kmers_default = -1;  // OFF affects processing
  static constexpr auto do_length_default = 1;
  static constexpr auto do_n_content_default = 1;
  static constexpr auto do_overrep_default = 1;  // affects processing
  static constexpr auto do_qual_base_default = 1;
  static constexpr auto do_qual_seq_default = 1;
  static constexpr auto do_sequence_default = 1;
  static constexpr auto do_tiles_default = -1;  // OFF affects processing

  int do_groups_{};

  int do_adap_{};
  int do_dups_{};
  int do_gc_content_{};
  int do_kmers_{};
  int do_length_{};
  int do_n_content_{};
  int do_overrep_{};
  int do_qual_base_{};
  int do_qual_seq_{};
  int do_sequence_{};
  int do_tiles_{};

  static std::vector<std::string> labels;
};

#endif  // SRC_RUN_MODE_HPP_
