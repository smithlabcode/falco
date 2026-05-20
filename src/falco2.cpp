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
#include "bam_file.hpp"
#include "duplication_results.hpp"
#include "falco_file_format.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "format_output.hpp"
#include "kmer_counter.hpp"
#include "quality_score.hpp"
#include "tile_processor.hpp"

#include "CLI11/CLI11.hpp"

#include <config.h>

#include <algorithm>
#include <array>
#include <compare>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <mutex>
#include <numeric>
#include <print>
#include <queue>
#include <ranges>
#include <stdexcept>
#include <stop_token>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

struct falco_results {
  static constexpr auto alphabet_size = 4;
  std::uint64_t n_reads{};
  std::uint64_t max_read_len{};
  std::vector<std::array<std::uint64_t, alphabet_size>> nucs;
  std::array<std::uint64_t, 101> gcs{};  // NOLINT (*-avoid-magic-numbers)
  std::vector<std::uint64_t> n_counts;
  std::vector<std::uint64_t> lengths;
  std::vector<std::array<std::uint64_t, falco::max_qual_val>> qual_by_pos;
  std::array<std::uint64_t, falco::max_qual_val> qual_by_read{};
  duplication_results dr;
  adapter_matcher am;
  std::string seq;
  std::string filename;

  falco_results() : lengths(1, 0) {}  // in case all reads have length 0

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
  process_one_read(this self_t &self, const auto &rec) -> void {
    self.process_one_read_impl(rec);
  }

  template <typename rec_t, typename self_t>
  auto
  process_reads(this self_t &self, auto cursor, const auto lim) {
    rec_t rec{};
    while (cursor < lim && (rec = get_next(cursor, lim))) {
      self.process_one_read(rec);
      ++self.n_reads;
    }
  }

  // clang-format off
  [[nodiscard]] auto make_seq_begin(const fqrec &rec) { return get_seq(rec); }
  [[nodiscard]] auto get_seq_begin(const fqrec &rec) { return get_seq(rec); }
  // clang-format on

  [[nodiscard]] auto
  make_seq_begin(const bamrec &rec) {
    auto itr = std::begin(seq);
    auto rec_seq_itr = get_seq(rec);
    const auto rec_seq_end = get_seq_end(rec);
    while (rec_seq_itr != rec_seq_end)
      *itr++ = *rec_seq_itr++;
    if (rec.is_rev)
      // ADS: only revcomp the prefix that will be used for this read
      revcomp(std::begin(seq), std::begin(seq) + get_seq_size(rec));
    return std::data(seq);
  }

  [[nodiscard]] auto
  get_seq_begin([[maybe_unused]] const bamrec &rec) {
    return std::cbegin(seq);
  }

  [[nodiscard]] auto
  process_qual(const bamrec &rec) {
    return rec.is_rev
             ? count_quals_rev(get_qual(rec), get_qual_end(rec), qual_by_pos)
             : count_quals(get_qual(rec), get_qual_end(rec), qual_by_pos);
  }

  [[nodiscard]] auto
  process_qual(const fqrec &rec) {
    return count_quals(get_qual(rec), get_qual_end(rec), qual_by_pos);
  }

  auto
  process_one_read_impl(const auto &rec) {
    // NOLINTBEGIN (*-pro-bounds-constant-array-index)
    static constexpr auto pct_int = [](const auto a, const auto b) {
      return (100 * a) / b;  // NOLINT (cppcoreguidelines-avoid-magic-numbers)
    };
    const auto read_len = static_cast<std::uint32_t>(get_seq_size(rec));
    if (read_len > max_read_len) {
      resize(read_len);
      if constexpr (std::is_same_v<std::decay_t<decltype(rec)>, bamrec>)
        seq.resize(read_len);
    }
    ++lengths[read_len];
    if (read_len == 0) [[unlikely]]
      return;
    max_read_len = read_len > max_read_len ? read_len : max_read_len;
    const auto seq_itr = make_seq_begin(rec);
    const auto seq_end = seq_itr + read_len;
    count_nucs(seq_itr, seq_end, nucs);
    const auto gc = count_gc(seq_itr, seq_end);
    ++gcs[pct_int(gc, read_len)];
    count_ns(seq_itr, seq_end, n_counts);
    const auto tot = process_qual(rec);
    ++qual_by_read[tot / read_len];
    count_seqs(seq_itr, read_len, dr);
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
    auto nucs_no_n = nucs;
    for (auto i = 0u; i < std::size(nucs_no_n); ++i)
      nucs_no_n[i][3] -= n_counts[i];
    return nucs_no_n;
  };

  auto
  fix_bam_qual_encoding() {  // BAM has no +33 for printable encoding
    // cppcheck-suppress constParameterReference
    const auto shift_and_fill = [&](auto &x) {
      std::shift_right(std::begin(x), std::end(x), falco::bam_qual_offset);
      std::ranges::fill_n(std::begin(x), falco::bam_qual_offset, 0);
    };
    shift_and_fill(qual_by_read);
    std::ranges::for_each(qual_by_pos, shift_and_fill);
  };

  template <typename self_t>
  [[nodiscard]] auto
  string(this self_t &self) {
    return self.string_impl();
  }

  [[nodiscard]] auto
  string_impl() const {
    const auto nucs_no_n = fix_nucs_for_ns();
    const auto total_nucs = tabular_dot(lengths);
    const auto gc_acc = [](const auto a, const auto &nuc) {
      return a + nuc[1] + nuc[3];
    };
    const auto total_gc = std::accumulate(std::cbegin(nucs_no_n),
                                          std::cend(nucs_no_n), 0ul, gc_acc);
    const auto encoding = identify_quality_score_encoding(qual_by_pos);
    const auto gt0 = [](const auto c) { return c > 0; };
    const std::uint64_t min_read_len =
      std::distance(std::cbegin(lengths), std::ranges::find_if(lengths, gt0));
    auto r = format_basic_stats(filename, n_reads, min_read_len, max_read_len,
                                total_gc, total_nucs, encoding);
    const auto qual_offset = identify_quality_score_offset(qual_by_pos);
    r += format_qual_by_pos(qual_by_pos, qual_offset);
    r += format_qual_by_read(qual_by_read, qual_offset);
    r += format_base_composition(nucs_no_n);     // base composition
    r += format_gc_content(gcs);                 // GC content
    r += format_n_counts(n_counts, nucs);        // N content
    r += format_read_lengths(lengths);           // read lengths
    r += dr.format_duplication_levels(n_reads);  // duplication results
    r += dr.format_overrepresented(n_reads);     // overrepresented sequences
    r += am.string(n_reads, max_read_len);       // adapter content
    return r;
  }
};

struct falco_results_tile : public falco_results {
  tile_processor tp;

  auto
  process_one_read_impl(const auto &rec) {
    falco_results::process_one_read_impl(rec);
    if (n_reads == tp.next_tile_read) {
      if (max_read_len > tp.max_read_len)
        tp.resize(max_read_len);
      tp.update_tile_id(get_name(rec), get_name_end(rec));
      tp(rec);
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
  string_impl() const {
    return falco_results::string_impl() + tp.string(max_read_len);
  }
};

struct falco_results_kmer : public falco_results {
  kmer_counter kc;

  auto
  process_one_read_impl(const auto &rec) {
    falco_results::process_one_read_impl(rec);
    if (n_reads == kc.next_kmer_read) {
      if (max_read_len > kc.max_read_len)
        kc.resize(max_read_len);
      kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
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
  string_impl() const {
    return falco_results::string_impl() + kc.string(n_reads);
  }
};

struct falco_results_tile_kmer : public falco_results_tile {
  kmer_counter kc;

  auto
  process_one_read_impl(const auto &rec) {
    falco_results_tile::process_one_read_impl(rec);
    if (n_reads == kc.next_kmer_read) {
      if (max_read_len > kc.max_read_len)
        kc.resize(max_read_len);
      kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
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
  string_impl() const {
    return falco_results_tile::string_impl() + kc.string(n_reads);
  }
};

template <typename results_t, typename rec_t> struct thread_pool {
  using task_t = std::pair<typename rec_t::pos_t, typename rec_t::pos_t>;
  std::mutex task_available_mtx;
  std::condition_variable_any task_available;
  std::mutex n_tasks_mtx;
  std::condition_variable finished;
  std::uint32_t n_tasks{};  // unfinished tasks (not same as unstarted)
  std::queue<task_t> tasks;
  std::vector<results_t> results;
  std::vector<std::jthread> workers;

  explicit thread_pool(std::uint32_t n_threads) : results(n_threads) {
    workers.reserve(n_threads);
    for (auto th_id = 0u; th_id < n_threads; ++th_id) {
      workers.emplace_back([th_id, this](const std::stop_token &stop) {
        auto &r = results[th_id];
        while (true) {
          std::pair<typename rec_t::pos_t, typename rec_t::pos_t> task;
          {
            std::unique_lock l(task_available_mtx);
            task_available.wait(l, stop, [this] { return !tasks.empty(); });
            if (stop.stop_requested() && tasks.empty())
              return;
            task = std::move(tasks.front());
            tasks.pop();
          }
          r.template process_reads<rec_t>(task.first, task.second);
          if (std::lock_guard l(n_tasks_mtx); --n_tasks == 0)
            finished.notify_all();
        }
      });
    }
  }

  void
  push_task(const task_t &chunk) {
    {
      std::scoped_lock lk(task_available_mtx, n_tasks_mtx);
      tasks.emplace(chunk);
      ++n_tasks;
    }
    task_available.notify_one();
  }

  void
  wait() {
    std::unique_lock l(n_tasks_mtx);
    finished.wait(l, [this] { return n_tasks == 0; });
  }
};

template <typename results_t>
static inline auto
accumulate_results(std::vector<results_t> &r) {
  // pointer jumping strategy
  auto id = std::views::iota(0, std::ssize(r)) | std::ranges::to<std::vector>();
  while (std::size(id) > 1) {
    {
      std::vector<std::jthread> workers;
      for (auto i = 0u; i + 1 < std::size(id); i += 2)
        // NOLINTNEXTLINE (performance-inefficient-vector-operation)
        workers.emplace_back([&, i] { r[id[i]] += r[id[i + 1]]; });
    }
    auto j = 0;
    for (auto i = 0; i + 1 < std::ssize(id); i += 2)
      id[j++] = id[i];
    if (std::ssize(id) % 2 == 1)
      id[j++] = id.back();
    id.resize(j);
  }
}

template <typename results_t>
static auto
run(const std::string &infile, auto &reads_file, const auto n_threads,
    const auto &outfile) {
  using rec_t = std::decay_t<decltype(reads_file)>::rec_t;
  thread_pool<results_t, rec_t> tpool(n_threads);
  while (reads_file) {
    reads_file.load_next();
    std::ranges::for_each(get_chunks(reads_file, n_threads),
                          [&](const auto &c) { tpool.push_task(c); });
    tpool.wait();
    // {
    //   std::vector<std::jthread> workers;
    //   for (auto th_id = 0u; th_id < n_threads; ++th_id)
    //     workers.emplace_back([&, th_id] {
    //       fr[th_id].template process_reads<rec_t>(chunks[th_id].first,
    //                                               chunks[th_id].second);
    //     });
    // }
  }
  std::ofstream out(outfile);
  if (!out)
    throw std::runtime_error("failed to open file: " + outfile);

  accumulate_results(tpool.results);
  auto &results = tpool.results.front();
  if constexpr (std::is_same_v<std::decay_t<decltype(reads_file)>, bam_file>)
    results.fix_bam_qual_encoding();
  results.filename = infile;
  std::print(out, "{}", results.string());
}

static auto
run_selector(const run_mode mode, const std::string &infile, auto &reads_file,
             const auto n_threads, const auto &outfile) {
  if (tiles(mode) && kmers(mode))
    run<falco_results_tile_kmer>(infile, reads_file, n_threads, outfile);
  else if (tiles(mode))
    run<falco_results_tile>(infile, reads_file, n_threads, outfile);
  else if (kmers(mode))
    run<falco_results_kmer>(infile, reads_file, n_threads, outfile);
  else
    run<falco_results>(infile, reads_file, n_threads, outfile);
}

int
main(int argc, char *argv[]) {
  try {
    static constexpr auto buf_size_defulat = 512 * 1024 * 1024;
    std::string infile;
    std::string outfile;
    std::int32_t buf_size{buf_size_defulat};
    std::uint32_t n_threads{1};

    bool do_tiles{};
    bool do_kmers{};
    bool verbose{};

    CLI::App app{PROJECT_NAME};
    argv = app.ensure_utf8(argv);
    if (argc >= 2)
      app.footer(PROJECT_NAME);

    // clang-format off
    app.set_help_flag("-h,--help", "Print a detailed help message and exit");
    app.add_option("-i,--input", infile, "FASTQ filename")
      ->required()
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "output filename")
      ->required();
    app.add_option("-s,--size", buf_size, "buffer size");
    app.add_option("-t,--threads", n_threads, "number of threads");
    app.add_flag("-v,--verbose", verbose, "print more run info");
    app.add_flag("--tiles,!--no-tiles", do_tiles,
                 std::format("toggle analysis per tiles (default: {})", do_tiles));
    app.add_flag("--kmers,!--no-kmers", do_kmers,
                 std::format("toggle analysis for kmers (default: {})", do_kmers));
    // clang-format on

    if (argc < 2) {
      std::cout << app.help() << '\n';
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const auto [infmt, infmt_descr] = get_file_format(infile);
    if (verbose)
      std::println("input file: {}\n"
                   "input file format: {}",
                   infile, infmt_descr);

    const bool has_tiles = tile_processor::set_preceding_colons(infile);
    if (verbose && do_tiles && !has_tiles)
      std::println("running without tiles: tile ids not found in data");
    do_tiles = do_tiles && has_tiles;

    run_mode mode;
    mode.tiles(do_tiles);
    mode.kmers(do_kmers);

    if (infmt == file_format::bam) {
      bam_file reads_file(infile, buf_size, n_threads);
      run_selector(mode, infile, reads_file, n_threads, outfile);
    }
    else if (infmt == file_format::fastq_gz) {
      fastq_gz_file reads_file(infile, buf_size, n_threads);
      run_selector(mode, infile, reads_file, n_threads, outfile);
    }
    else if (infmt == file_format::fastq) {
      fastq_file reads_file(infile, buf_size);
      run_selector(mode, infile, reads_file, n_threads, outfile);
    }
    else {
      std::println("unsupported file format: {}", infmt_descr);
    }
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
