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

#include "adapter_matcher.hpp"
#include "duplication_results.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "fastq_record.hpp"
#include "kmer_counter.hpp"
#include "tile_processor.hpp"

#include "CLI11/CLI11.hpp"

#include <config.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <format>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <numeric>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <vector>

struct falco_results {
  static constexpr auto init_read_len = 256;
  static constexpr auto alphabet_size = 4;

  std::uint64_t n_reads{};
  std::uint64_t max_read_len{};
  std::uint64_t gc{};
  std::vector<std::array<std::uint64_t, alphabet_size>> nucs;
  std::array<std::uint64_t, 101> gcs{};  // NOLINT (*-avoid-magic-numbers)
  std::vector<std::uint64_t> n_counts;
  std::vector<std::uint64_t> lengths;
  std::vector<std::array<std::uint64_t, falco::max_qual_val>> qual_by_pos;
  std::array<std::uint64_t, falco::max_qual_val> qual_by_read{};
  duplication_results dr;
  adapter_matcher am;

  // clang-format off
  falco_results() :
    nucs(init_read_len, std::array<std::uint64_t, alphabet_size>{}),
    n_counts(init_read_len, 0),
    lengths(init_read_len + 1, 0),
    qual_by_pos(init_read_len, std::array<std::uint64_t, falco::max_qual_val>{}),
    am(init_read_len)
  {}
  // clang-format on

  auto
  resize(const std::uint32_t updated_length) {
    nucs.resize(updated_length);
    n_counts.resize(updated_length);
    lengths.resize(updated_length + 1);  // need one extra here
    qual_by_pos.resize(updated_length);
    am.resize(updated_length);
  }

  template <typename self_t>
  auto
  process_one_read(this self_t &&self, const auto &reads_buf,
                   const auto &rec) -> void {
    self.process_one_read_impl(reads_buf, rec);
  }

  template <typename self_t>
  auto
  process_reads(this self_t &&self, const auto &reads_buf, std::int64_t cursor,
                const std::int64_t lim) {
    using rec_t = std::decay_t<decltype(reads_buf)>::rec_t;
    rec_t rec{};
    while (cursor < lim && (rec = get_next(reads_buf, cursor, lim))) {
      self.process_one_read(reads_buf, rec);
      ++self.n_reads;
    }
  }

  auto
  process_one_read_impl(const auto &reads_buf, const auto &rec) {
    // NOLINTBEGIN (*-pro-bounds-constant-array-index)
    static constexpr auto pct_int = [](const auto a, const auto b) {
      return (100 * a) / b;  // NOLINT (cppcoreguidelines-avoid-magic-numbers)
    };
    const auto read_len = get_seq_size(rec);
    const auto seq_itr = get_seq(reads_buf, rec);
    const auto seq_end = get_seq_end(reads_buf, rec);
    if (read_len > max_read_len)
      resize(read_len);
    max_read_len = read_len > max_read_len ? read_len : max_read_len;
    ++lengths[read_len];
    count_nucs(seq_itr, seq_end, nucs);
    const auto curr_gc = count_gc(seq_itr, seq_end);
    ++gcs[pct_int(curr_gc, read_len)];
    gc += curr_gc;
    count_ns(seq_itr, seq_end, n_counts);
    const auto qual_itr = get_qual(reads_buf, rec);
    const auto qual_end = get_qual_end(reads_buf, rec);
    const auto qtot = count_quals(qual_itr, qual_end, qual_by_pos) / read_len;
    ++qual_by_read[qtot];
    count_seqs(seq_itr, read_len, n_reads, dr);
    am.match_adapters(seq_itr, read_len);
    // NOLINTEND (*-pro-bounds-constant-array-index)
  }

  auto
  operator+=(const falco_results &rhs) -> const falco_results & {
    n_reads += rhs.n_reads;
    max_read_len = std::max(max_read_len, rhs.max_read_len);
    two_dim_add(nucs, rhs.nucs);
    add(gcs, rhs.gcs);
    vec_add(lengths, rhs.lengths);
    vec_add(n_counts, rhs.n_counts);
    two_dim_add(qual_by_pos, rhs.qual_by_pos);
    add(qual_by_read, rhs.qual_by_read);

    dr += rhs.dr;
    am += rhs.am;

    return *this;
  }

  [[nodiscard]] auto
  fix_nucs_for_ns() const {  // Ns were counted among the C in nucs
    auto x = nucs;
    for (auto i = 0; i < std::size(x); ++i)
      x[i][3] -= n_counts[i];
    return x;
  };

  template <typename self_t>
  [[nodiscard]] auto
  string(this const self_t &self) -> std::string {
    return self.string_impl();
  }

  [[nodiscard]] auto
  string_impl() const -> std::string {
    const auto nucs_no_n = fix_nucs_for_ns();
    const auto total_nucs = tabular_dot(lengths);
    const auto total_gc = std::accumulate(
      std::cbegin(nucs_no_n), std::cend(nucs_no_n), 0ul,
      [](const auto a, const auto &nuc) { return a + nuc[1] + nuc[3]; });

    auto r = format_basic_stats(/*filename, */ n_reads, max_read_len, total_gc,
                                total_nucs);
    r += format_qual_by_pos(qual_by_pos, max_read_len);     // qual by pos
    r += format_qual_by_read(qual_by_read, max_read_len);   // qual by read
    r += format_base_composition(nucs_no_n, max_read_len);  // base composition
    r += format_gc_content(gcs, max_read_len);              // GC content
    r += format_n_counts(n_counts, nucs, max_read_len);     // N content
    r += format_read_lengths(lengths, max_read_len);        // read lengths

    r += dr.string();                         // duplication results
    r += dr.format_overrepresented(n_reads);  // overrepresented sequences
    r += am.string(n_reads, max_read_len);    // adapter content

    return r;
  }
};

struct falco_results_tile : public falco_results {
  tile_processor tp;

  falco_results_tile() : tp(init_read_len) {}

  auto
  resize(const std::uint32_t updated_length) {
    falco_results::resize(updated_length);
    tp.resize(updated_length);
  }

  auto
  process_one_read_impl(const auto &reads_buf, const auto &rec) {
    falco_results::process_one_read_impl(reads_buf, rec);
    if (n_reads == tp.next_tile_read) {
      const auto name_itr = get_name(reads_buf, rec);
      const auto name_end = get_name_end(reads_buf, rec);
      tp.update_tile_id(name_itr, name_end);
      const auto qual_itr = get_qual(reads_buf, rec);
      const auto qual_end = get_qual_end(reads_buf, rec);
      tp(qual_itr, qual_end);
      tp.next_tile_read += tp.tile_step;
    }
  }

  auto
  operator+=(const falco_results_tile &rhs) -> const falco_results_tile & {
    falco_results::operator+=(rhs);
    tp += rhs.tp;
    return *this;
  }

  [[nodiscard]] auto
  string_impl() const -> std::string {
    return falco_results::string_impl() + tp.string(max_read_len);
  }
};

struct falco_results_kmer : public falco_results {
  kmer_counter kc;

  falco_results_kmer() : kc(init_read_len) {}

  auto
  resize(const std::uint32_t updated_length) {
    falco_results::resize(updated_length);
    kc.resize(updated_length);
  }

  auto
  process_one_read_impl(const auto &reads_buf, const auto &rec) {
    falco_results::process_one_read_impl(reads_buf, rec);
    if (n_reads == kc.next_kmer_read) {
      const auto read_len = get_seq_size(rec);
      const auto seq = get_seq(reads_buf, rec);
      kc.count_kmers(n_reads, seq, read_len);
      kc.next_kmer_read += kmer_counter::kmer_step;
    }
  }

  auto
  operator+=(const falco_results_kmer &rhs) -> const falco_results_kmer & {
    falco_results::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  string_impl() const -> std::string {
    return falco_results::string_impl() + kc.string(n_reads);
  }
};

struct falco_results_tile_kmer : public falco_results_tile {
  kmer_counter kc;

  falco_results_tile_kmer() : kc(init_read_len) {}

  auto
  resize(const std::uint32_t updated_length) {
    falco_results_tile::resize(updated_length);
    kc.resize(updated_length);
  }

  auto
  process_one_read_impl(const auto &reads_buf, const auto &rec) {
    falco_results_tile::process_one_read_impl(reads_buf, rec);
    if (n_reads == kc.next_kmer_read) {
      const auto read_len = get_seq_size(rec);
      const auto seq = get_seq(reads_buf, rec);
      kc.count_kmers(n_reads, seq, read_len);
      kc.next_kmer_read += kmer_counter::kmer_step;
    }
  }

  auto
  operator+=(const falco_results_tile_kmer &rhs)
    -> const falco_results_tile_kmer & {
    falco_results_tile::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  string_impl() const -> std::string {
    return falco_results_tile::string_impl() + kc.string(n_reads);
  }
};

template <typename results_t>
static auto
run(auto &reads_file, const auto n_threads, const auto &output_filename) {
  std::vector<results_t> fr(n_threads);
  while (reads_file) {
    auto reads_buf = reads_file.get_next();
    const auto chunks = get_chunks(reads_file, reads_buf, n_threads);
    {
      std::vector<std::jthread> workers;
      for (auto th_id = 0; th_id < n_threads; ++th_id)
        // NOLINTNEXTLINE (performance-inefficient-vector-operation)
        workers.emplace_back([&, th_id] {
          fr[th_id].process_reads(reads_buf, chunks[th_id].first,
                                  chunks[th_id].second);
        });
    }
  }
  std::ofstream out(output_filename);
  if (!out)
    throw std::runtime_error("failed to open file: " + output_filename);
  const auto results = std::reduce(std::cbegin(fr), std::cend(fr));
  std::println(out, "{}", results.string());
}

int
main(int argc, char *argv[]) {
  try {
    static constexpr auto buf_size_defulat = 512 * 1024 * 1024;
    std::string fastq_filename;
    std::string output_filename;
    std::int64_t buf_size{buf_size_defulat};
    std::int64_t n_threads{1};

    bool do_tiles{};
    bool do_kmers{};
    bool is_gzip{};

    CLI::App app{PROJECT_NAME};
    argv = app.ensure_utf8(argv);
    if (argc >= 2)
      app.footer(PROJECT_NAME);

    // clang-format off
    app.set_help_flag("-h,--help", "Print a detailed help message and exit");
    app.add_option("-i,--input", fastq_filename, "FASTQ filename")
      ->required()
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", output_filename, "output filename")
      ->required();
    app.add_option("-s,--size", buf_size, "buffer size");
    app.add_option("-t,--threads", n_threads, "number of threads");
    app.add_flag("--tiles", do_tiles, "report results per tile");
    app.add_flag("--kmers", do_kmers, "report results for kmers");
    app.add_flag("-z,--zip", is_gzip, "input file is gzip");
    // clang-format on

    if (argc < 2) {
      std::cout << app.help() << '\n';
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const bool has_tiles = tile_processor::set_preceding_colons(fastq_filename);
    do_tiles = do_tiles && has_tiles;

    if (is_gzip) {
      fastq_gz_file reads_file(fastq_filename, buf_size);
      if (do_tiles && do_kmers)
        run<falco_results_tile_kmer>(reads_file, n_threads, output_filename);
      else if (do_tiles)
        run<falco_results_tile>(reads_file, n_threads, output_filename);
      else if (do_kmers)
        run<falco_results_kmer>(reads_file, n_threads, output_filename);
      else
        run<falco_results>(reads_file, n_threads, output_filename);
    }
    else {
      fastq_file reads_file(fastq_filename, buf_size);
      if (do_tiles && do_kmers)
        run<falco_results_tile_kmer>(reads_file, n_threads, output_filename);
      else if (do_tiles)
        run<falco_results_tile>(reads_file, n_threads, output_filename);
      else if (do_kmers)
        run<falco_results_kmer>(reads_file, n_threads, output_filename);
      else
        run<falco_results>(reads_file, n_threads, output_filename);
    }
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
