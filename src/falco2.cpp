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

// clang-format off
static constexpr auto about =
  R"(Falco2: an in-progress redesign and rewrite of falco)";

static constexpr auto description =
  R"(Falco2 aims to improve speed while still doing the same analysis as in falco.
There will likely be changes to some of the statistics, including the way
read duplication is analyzed (borrowing from preseq).
)";
// clang-format on

#include "adapter_matcher.hpp"
#include "bam_file.hpp"
#include "contaminants.hpp"
#include "duplication_results.hpp"
#include "falco_file_format.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "format_output.hpp"
#include "kmer_counter.hpp"
#include "quality_score.hpp"
#include "tile_processor.hpp"

#include "CLI11/CLI11.hpp"
#include "nlohmann/json.hpp"

#include <config.h>

#include <license.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <compare>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <format>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <mutex>
#include <new>
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

struct file_info {
  std::string name;
  falco::file_format format{};
  std::string description;
  std::uint64_t size{};
  falco::encoding encoding{};
  bool has_tiles{};

  [[nodiscard]] auto
  string() const -> std::string {
    static constexpr auto n_indent = 4;
    nlohmann::json data = *this;
    return data.dump(n_indent);
  }

  NLOHMANN_DEFINE_TYPE_INTRUSIVE(file_info, name, format, description, size,
                                 encoding, has_tiles);
};

struct alignas(std::hardware_destructive_interference_size) falco_results {
  std::uint64_t n_reads{};
  std::uint64_t max_read_len{};
  std::vector<falco::nuc_array> nucs;
  falco::gc_content_array gcs{};
  std::vector<std::uint64_t> n_counts;
  std::vector<std::uint64_t> lengths;
  std::vector<falco::qual_array> qual_by_pos;
  falco::qual_array qual_by_read{};
  duplication_results dr;
  adapter_matcher am;
  std::string seq;

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
  finalize_qual_encoding(this self_t &self, const auto enc) {
    self.finalize_qual_encoding_impl(enc);
  }

  auto
  finalize_qual_encoding_impl(const auto enc) {
    adjust_fastq_qual_encoding(qual_by_pos, qual_by_read, enc);
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
  [[nodiscard]] auto get_seq_begin(const fqrec &rec) { return get_seq(rec); }
  [[nodiscard]] auto make_seq_begin(const fqrec &rec) { return get_seq(rec); }
  [[nodiscard]] auto get_seq_begin(const bamrec &) { return std::data(seq); }
  // clang-format on

  [[nodiscard]] auto
  make_seq_begin(const bamrec &rec) {
    const auto complement = [](const auto a) {
      return "TNGNNNCNNNNNNNNNNNNA"[a - 'A'];
    };
    auto rec_seq_itr = get_seq(rec);
    const auto rec_seq_end = get_seq_end(rec);
    if (rec.is_rev) {
      auto itr = std::begin(seq) + get_seq_size(rec);
      while (rec_seq_itr != rec_seq_end)
        *(--itr) = complement(*rec_seq_itr++);
    }
    else {
      auto itr = std::begin(seq);
      while (rec_seq_itr != rec_seq_end)
        *itr++ = *rec_seq_itr++;
    }
    return std::data(seq);
  }

  [[nodiscard]] auto
  process_quality_scores(const bamrec &rec) {
    return rec.is_rev
             ? count_quals_rev(get_qual(rec), get_qual_end(rec), qual_by_pos)
             : count_quals(get_qual(rec), get_qual_end(rec), qual_by_pos);
  }

  [[nodiscard]] auto
  process_quality_scores(const fqrec &rec) {
    return count_quals(get_qual(rec), get_qual_end(rec), qual_by_pos);
  }

  auto
  process_one_read_impl(const auto &rec) {
    // NOLINTBEGIN (cppcoreguidelines-pro-bounds-constant-array-index)
    static constexpr auto discrete_pct = [](const auto a, const auto b) {
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
    ++gcs[discrete_pct(gc, read_len)];
    count_ns(seq_itr, seq_end, n_counts);
    const auto tot = process_quality_scores(rec);
    ++qual_by_read[tot / read_len];
    dr.count_seqs(seq_itr, read_len);
    am.match_adapters(seq_itr, read_len);
    // NOLINTEND (cppcoreguidelines-pro-bounds-constant-array-index)
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
  adjust_nucs_for_ns() const {  // Ns were counted among the C in nucs
    auto nucs_no_n = nucs;
    for (auto i = 0u; i < std::size(nucs_no_n); ++i)
      nucs_no_n[i][3] -= n_counts[i];
    return nucs_no_n;
  };

  template <typename self_t>
  [[nodiscard]] auto
  string(this self_t &self, const file_info &info) {
    return self.string_impl(info);
  }

  [[nodiscard]] auto
  string_impl(const file_info &info) const {
    const auto nucs_no_n = adjust_nucs_for_ns();
    const auto total_nucs = tabular_dot(lengths);
    const auto gc_acc = [](const auto a, const auto &nuc) {
      return a + nuc[1] + nuc[3];  // NOLINT (*-avoid-magic-numbers)
    };
    const auto total_gc = std::accumulate(std::cbegin(nucs_no_n),
                                          std::cend(nucs_no_n), 0ul, gc_acc);
    const auto gt0 = [](const auto c) { return c > 0; };
    const std::uint64_t min_read_len =
      std::distance(std::cbegin(lengths), std::ranges::find_if(lengths, gt0));
    const auto encoding_label = get_quality_score_label(info.encoding);
    auto r = format_basic_stats(info.name, n_reads, min_read_len, max_read_len,
                                total_gc, total_nucs, encoding_label);
    r += format_qual_by_pos(qual_by_pos);
    r += format_qual_by_read(qual_by_read);
    r += format_base_composition(nucs_no_n);  // base composition
    r += format_gc_content(gcs);              // GC content
    r += format_n_counts(n_counts, nucs);     // N content
    r += format_read_lengths(lengths);        // read lengths
    r += dr.format_duplication_levels();      // duplication results
    r += dr.format_overrepresented();         // overrepresented sequences
    r += am.string(n_reads);                  // adapter content
    return r;
  }
};

struct falco_results_tile : public falco_results {
  tile_processor tp;

  auto
  finalize_qual_encoding_impl(const auto enc) {
    falco_results::finalize_qual_encoding_impl(enc);
    tp.trim();
    tp.adjust_fastq_qual_encoding(enc);
  }

  auto
  process_one_read_impl(const auto &rec) {
    falco_results::process_one_read_impl(rec);
    tp(rec);
  }

  auto
  operator+=(const falco_results_tile &rhs) -> const falco_results_tile & {
    falco_results::operator+=(rhs);
    tp += rhs.tp;
    return *this;
  }

  [[nodiscard]] auto
  string_impl(const file_info &info) const {
    return falco_results::string_impl(info) + tp.string(max_read_len);
  }
};

struct falco_results_kmer : public falco_results {
  kmer_counter kc;

  auto
  process_one_read_impl(const auto &rec) {
    falco_results::process_one_read_impl(rec);
    kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
  }

  auto
  operator+=(const falco_results_kmer &rhs) -> const falco_results_kmer & {
    falco_results::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  string_impl(const file_info &info) const {
    return falco_results::string_impl(info) + kc.string();
  }
};

struct falco_results_tile_kmer : public falco_results_tile {
  kmer_counter kc;

  auto
  finalize_qual_encoding_impl(const auto enc) {
    falco_results_tile::finalize_qual_encoding_impl(enc);
  }

  auto
  process_one_read_impl(const auto &rec) {
    falco_results_tile::process_one_read_impl(rec);
    kc.count_kmers(get_seq_begin(rec), get_seq_size(rec));
  }

  auto
  operator+=(const falco_results_tile_kmer &rhs)
    -> const falco_results_tile_kmer & {
    falco_results_tile::operator+=(rhs);
    kc += rhs.kc;
    return *this;
  }

  [[nodiscard]] auto
  string_impl(const file_info &info) const {
    return falco_results_tile::string_impl(info) + kc.string();
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
run(file_info &info, auto &reads_file, const auto n_threads,
    const auto &outfile) {
  static constexpr auto n_chunks_per_thread = 2;
  using rec_t = std::decay_t<decltype(reads_file)>::rec_t;
  thread_pool<results_t, rec_t> tpool(n_threads);
  while (reads_file) {
    reads_file.load_next();
    std::ranges::for_each(
      get_chunks(reads_file, n_threads * n_chunks_per_thread),
      [&](const auto &c) { tpool.push_task(c); });
    tpool.wait();
    // ADS: below is how we'd do it without a pool, and it's likely as fast with
    // current bottlencks being elsewhere.

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

  set_quality_score_encoding(results.qual_by_pos, info);
  if (!is_mapped_reads(info.format))
    results.finalize_qual_encoding(info.encoding);
  std::print(out, "{}", results.string(info));
}

static auto
run_mode_selector(const run_mode mode, file_info &info, auto &reads_file,
                  const auto n_threads, const auto &outfile) {
  if (tiles(mode) && kmers(mode))
    run<falco_results_tile_kmer>(info, reads_file, n_threads, outfile);
  else if (tiles(mode))
    run<falco_results_tile>(info, reads_file, n_threads, outfile);
  else if (kmers(mode))
    run<falco_results_kmer>(info, reads_file, n_threads, outfile);
  else
    run<falco_results>(info, reads_file, n_threads, outfile);
}

int
main(int argc, char *argv[]) {
  try {
    static constexpr auto buf_size_default = 256 * 1024 * 1024;
    std::string infile;
    std::string contam_file;
    std::string outfile;
    std::int64_t buf_size{buf_size_default};
    std::uint32_t n_threads{1};

    bool do_tiles{};
    bool do_kmers{};
    bool verbose{};

    using std::literals::string_literals::operator""s;
    const auto size_from_units = CLI::AsNumberWithUnit(std::map{
      std::pair{"G"s, gigabytes},
      {"M"s, megabytes},
      {"k"s, kilobytes},
    });

    CLI::App app{about};
    argv = app.ensure_utf8(argv);
    app.usage(
      std::format("Usage: {} [options] -o OUTFILE -i INFILE", PROJECT_NAME));
    if (argc >= 2)
      app.footer(description);

    // NOLINTNEXTLINE (cppcoreguidelines-avoid-magic-numbers)
    app.get_formatter()->long_option_alignment_ratio(0.2);
    app.set_help_flag("-h,--help", "Print more detailed help");
    app.set_version_flag("--version", VERSION, "Print program version");
    // clang-format off
    app.add_flag("--license",
                 [&](auto) { std::print("{}", license_text); throw CLI::Success(); },
                 "Print full license")
      ->callback_priority(CLI::CallbackPriority::First);
    // clang-format on
    app
      .add_option("-i,--input", infile,
                  "Input file: FASTQ (plain, gz or bgzf) or BAM/SAM")
      ->required()
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outfile, "Output file")
      ->required()
      ->option_text("FILE");
    app
      .add_option("-c,--contaminants", contam_file,
                  "File of contaminant sequences to use")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app
      .add_option("-t,--threads", n_threads,
                  std::format("Threads to use (this machine supports: {})",
                              std::thread::hardware_concurrency()))
      ->option_text(std::format("[{}]", n_threads));
    app
      .add_option("-m,--mem", buf_size,
                  "Memory buffer size for IO (G/M/K units ok)")
      ->option_text(std::format("[{}]", size_to_units(buf_size_default)))
      ->capture_default_str()
      ->transform(size_from_units);
    app.add_flag("-v,--verbose", verbose, "Print more info while running.");
    app.add_flag("--tiles", do_tiles, "Enable per-tile analysis");
    app.add_flag("--kmers", do_kmers, "Enable k-mer analysis");

    const auto start_time{std::chrono::high_resolution_clock::now()};

    if (argc < 2) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    if (n_threads > std::thread::hardware_concurrency())
      n_threads = std::thread::hardware_concurrency();

    if (verbose)
      std::print("input file: {}\n"
                 "output file: {}\n"
                 "memory requested: {}\n"
                 "threads requested: {}\n"
                 "tile analysis requested: {}\n"
                 "k-mer analysis requested: {}\n",
                 infile, outfile, size_to_units(buf_size), n_threads, do_tiles,
                 do_kmers);

    if (!contam_file.empty()) {
      load_contaminants(contam_file);
      if (verbose)
        std::print("contaminants file: {}\n"
                   "number of contaminants: {}\n",
                   contam_file, std::size(contaminants));
    }

    const auto [input_format, format_description] = get_file_format(infile);
    if (verbose)
      std::println("input file: {}\n"
                   "input file format: {}",
                   infile, format_description);

    const bool has_tiles = tile_processor::set_preceding_colons(infile);
    if (verbose && do_tiles && !has_tiles)
      std::println("running without tiles: tile ids not found in data");
    do_tiles = do_tiles && has_tiles;

    run_mode mode;
    mode.tiles(do_tiles);
    mode.kmers(do_kmers);

    const auto [est_n_reads, filesize] = [&] {
      if (input_format == falco::file_format::bam)
        return estimate_n_reads_bam(infile);
      if (input_format == falco::file_format::fastq_gz)
        return estimate_n_reads_fastq_gz(infile);
      if (input_format == falco::file_format::fastq)
        return estimate_n_reads_fastq(infile);
      std::unreachable();
    }();

    file_info info{
      .name = infile,
      .format = input_format,
      .description = format_description,
      .size = filesize,
      .has_tiles = has_tiles,
    };

    duplication_results::initialize(est_n_reads);

    if (is_mapped_reads(input_format)) {
      bam_file reads_file(infile, buf_size, n_threads);
      run_mode_selector(mode, info, reads_file, n_threads, outfile);
    }
    else if (input_format == falco::file_format::fastq_gz) {
      fastq_gz_file reads_file(infile, buf_size, n_threads);
      run_mode_selector(mode, info, reads_file, n_threads, outfile);
    }
    else if (input_format == falco::file_format::fastq) {
      fastq_file reads_file(infile, buf_size);
      run_mode_selector(mode, info, reads_file, n_threads, outfile);
    }
    else {
      std::println("unsupported file format: {}", format_description);
    }

    if (verbose) {
      // ADS: using 'count()' because macos has locale issues formatting times
      const auto dur = [](const auto d) {
        return std::chrono::duration_cast<std::chrono::duration<double>>(d)
          .count();
      };
      const auto stop_time{std::chrono::high_resolution_clock::now()};
      std::print("total run time: {:.6g}s\n", dur(stop_time - start_time));
    }
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
