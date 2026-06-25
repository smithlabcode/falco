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

#include "adapter_set.hpp"
#include "bam_file.hpp"
#include "contaminants.hpp"
#include "falco_analyzer.hpp"
#include "falco_file_format.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "quality_score.hpp"
#include "results_collector.hpp"
#include "thread_pool_wrapper.hpp"
#include "tile_processor.hpp"

#include "CLI11/CLI11.hpp"

#include <config.h>
#include <license.h>

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <format>
#include <fstream>
#include <functional>
#include <iterator>
#include <map>
#include <memory>
#include <print>
#include <ranges>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

struct thread_counter {
  std::uint32_t workers{};
  std::uint32_t readers{};
  std::uint32_t decomp{};

  auto
  initialize(const auto input_format, const auto n_files) {
    const std::uint32_t max_threads = std::thread::hardware_concurrency();
    workers = std::min(workers, max_threads);
    readers = std::min(readers, max_threads);
    decomp = std::min(decomp, max_threads);
    if (readers == 0)
      readers = workers < n_files ? workers : n_files;
    if (is_bgzf(input_format) && decomp == 0)
      decomp = workers;
    if (!is_bgzf(input_format))
      // ADS: need a warning on this
      decomp = 0;
  }
};

static auto
write_file(const auto &filename, const auto &data) {
  std::ofstream out(filename);
  if (!out)
    throw std::runtime_error("failed to open file: " + filename);
  std::print(out, "{}", data);
}

template <typename results_t>
static auto
run(const run_mode &mode, std::vector<file_info> &infos, auto &reads_files,
    const auto n_threads, const auto &outdirs) {
  using rec_t = std::decay_t<decltype(reads_files)>::value_type::rec_t;
  static constexpr auto report_filename = "fastqc_data.txt";
  static constexpr auto html_filename = "fastqc_report.html";
  static constexpr auto summary_filename = "summary.txt";

  const auto n_files = std::size(reads_files);

  analyzer_t<results_t, rec_t> analyzer(n_threads.workers, n_threads.readers,
                                        n_files, reads_files, infos);

  for (const auto file_id : std::views::iota(0u, n_files))
    accumulate_results(analyzer.results, file_id);

  for (const auto [results, info, outdir] :
       std::views::zip(analyzer.results.front(), infos, outdirs)) {
    set_quality_score_encoding(results.qual_by_pos, info);
    results.finalize(mode, info);
    const auto grades = results.get_grades();
    const auto report = results.get_report(mode, info, grades);

    const auto report_path =
      (std::filesystem::path{outdir} / report_filename).string();
    write_file(report_path, report);

    const auto html = results.get_html(mode, info, grades);
    const auto html_path =
      (std::filesystem::path{outdir} / html_filename).string();
    write_file(html_path, html);

    const auto summary_path =
      (std::filesystem::path{outdir} / summary_filename).string();
    write_file(summary_path, grades.summary(info.name));
  }
}

static auto
run_mode_selector(const run_mode &mode, std::vector<file_info> &infos,
                  auto &reads_files, const auto n_threads,
                  const auto &outdirs) {
  if (do_tiles(mode) && do_kmers(mode))
    run<results_collector_tile_kmer>(mode, infos, reads_files, n_threads,
                                     outdirs);
  else if (do_tiles(mode))
    run<results_collector_tile>(mode, infos, reads_files, n_threads, outdirs);
  else if (do_kmers(mode))
    return run<results_collector_kmer>(mode, infos, reads_files, n_threads,
                                       outdirs);
  else
    run<results_collector>(mode, infos, reads_files, n_threads, outdirs);
}

static auto
start_analysis(const run_mode &mode, const auto buf_size, const auto n_threads,
               std::vector<file_info> &infos, const auto &infiles,
               const auto &outdirs) {
  // ADS: this function is bloated mostly to avoid allowing default construction
  // or move assignment in the fastq_file, fastq_gz_file, etc.
  const auto input_format = infos.front().format;
  const auto n_infiles = std::size(infiles);
  if (is_mapped_reads(input_format)) {
    falco_thread_pool p(n_threads.decomp);
    std::vector<bam_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size, p); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outdirs);
  }
  else if (input_format == falco::file_format::fastq_bgzf) {
    falco_thread_pool p(n_threads.decomp);
    std::vector<fastq_bgzf_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size, p); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outdirs);
  }
  else if (input_format == falco::file_format::fastq_gz) {
    std::vector<fastq_gz_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outdirs);
  }
  else if (input_format == falco::file_format::fastq) {
    std::vector<fastq_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outdirs);
  }
  else {
    const auto format_description = infos.front().description;
    std::println("unsupported file format: {}", format_description);
  }
}

[[nodiscard]] static auto
get_file_info(const auto &infiles) {
  std::vector<file_info> infos;
  for (const auto [file_id, infile] : std::views::enumerate(infiles)) {
    const auto [input_format, format_description] = get_file_format(infile);
    const auto tile_id_position = get_tile_info(infile);
    const bool has_tiles = (tile_id_position != 0);
    const auto [n_reads_est, read_len_est, filesize] = [&] {
      if (input_format == falco::file_format::bam)
        return estimate_n_reads_bam(infile);
      if (input_format == falco::file_format::fastq_bgzf)
        return estimate_n_reads_fastq_bgzf(infile);
      if (input_format == falco::file_format::fastq_gz)
        return estimate_n_reads_fastq_gz(infile);
      if (input_format == falco::file_format::fastq)
        return estimate_n_reads_fastq(infile);
      std::unreachable();
    }();
    infos.push_back({
      .name = infile,
      .format = input_format,
      .description = format_description,
      .size = filesize,
      .n_reads_est = n_reads_est,
      .read_len_est = read_len_est,
      .has_tiles = has_tiles,
      .tile_id_position = tile_id_position,
    });
  }
  return infos;
}

[[nodiscard]] static auto
make_outdirs(const auto &ins, const auto &outdir) -> std::vector<std::string> {
  namespace fs = std::filesystem;
  fs::create_directory(outdir);
  const auto compose_dirname = [&](const auto &fname) {
    const auto without_path = fs::path{fname}.filename();
    return (fs::path{outdir} / remove_extension(without_path)).string();
  };
  const auto dnames = ins | std::views::transform(compose_dirname);
  std::ranges::for_each(dnames, [](const auto &d) { fs::create_directory(d); });
  return dnames | std::ranges::to<std::vector>();
}

int
main(int argc, char *argv[]) {
  try {
    static constexpr auto buffer_size_default = 256 * 1024 * 1024;
    std::vector<std::string> infiles;
    std::string contam_file;
    std::string config_file;
    std::string adapters_file;
    std::string outdir;
    std::int64_t buffer_size{buffer_size_default};

    thread_counter n_threads{1, 0, 0};

    bool do_tiles{};
    bool do_kmers{};
    bool do_groups{};
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
      std::format("Usage: {} [options] -o OUTDIR INFILES", PROJECT_NAME));
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
    app.add_option("INFILES", infiles,
                   "Input file: FASTQ (plain, GZIP or BGZF) or BAM/SAM")
      ->required()
      ->option_text("FILES")
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outdir, "Output directory")
      ->required()
      ->option_text("DIR");
    app.add_option("--config", config_file,
                   "Configuration file")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-c,--contaminants", contam_file,
                   "File of contaminant sequences to use")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-a,--adapters", adapters_file, "File of adapters sequences to use")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-t,--threads", n_threads.workers,
                   std::format("Threads for analysis (this machine supports: {})",
                               std::thread::hardware_concurrency()))
      ->option_text(std::format("[{}]", n_threads.workers));
    app.add_option("-r,--readers", n_threads.readers,
                   "Threads for reading input (default: one per input file)");
    app.add_option("-d,--decomp", n_threads.decomp,
                   "Threads for BAM/BGZF decompression (default: analysis threads)");
    app.add_option("-m,--mem", buffer_size,
                   "Memory buffer size for IO (G/M/K units ok)")
      ->option_text(std::format("[{}]", size_to_units(buffer_size_default)))
      ->capture_default_str()
      ->transform(size_from_units);
    app.add_flag("-v,--verbose", verbose, "Print more info while running.");
    app.add_flag("--tiles", do_tiles, "Enable per-tile analysis");
    app.add_flag("--kmers", do_kmers, "Enable k-mer analysis");
    app.add_flag("--groups", do_groups, "Use base groups in output");
    // clang-format on

    const auto start_time{std::chrono::high_resolution_clock::now()};

    if (argc < 2) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    if (!config_file.empty())
      grader_set::instance(config_file);

    const auto outdirs = make_outdirs(infiles, outdir);

    if (!contam_file.empty()) {
      load_contaminants(contam_file);
      if (verbose)
        std::print("contaminants file: {}\n"
                   "number of contaminants: {}\n",
                   contam_file, std::size(contaminants));
    }

    if (!adapters_file.empty()) {
      const auto &as = adapter_set::instance(adapters_file);
      if (verbose)
        std::print("adapters file: {}\n"
                   "number of adapters: {}\n",
                   adapters_file, adapter_set::n_adapters());
      if (const auto [is_valid, message] = as.validate(); !is_valid) {
        std::println("{}", message);
        return EXIT_FAILURE;
      }
    }

    auto infos = get_file_info(infiles);

    const auto format_description = infos.front().description;
    const auto input_format = infos.front().format;
    const auto ok_fmt = [&](const auto &x) { return x.format == input_format; };
    if (!std::ranges::all_of(infos, ok_fmt)) {
      std::println("unequal input formats [expected {}]", format_description);
      return EXIT_FAILURE;
    }

    n_threads.initialize(input_format, std::size(infiles));

    // restrict buffer size to avoid using a possibly harmful amount of memory
    const auto get_sz = [](const auto &i) { return i.size; };
    const auto max_sz = std::ranges::max(std::views::transform(infos, get_sz));
    buffer_size = buffer_size < max_sz ? buffer_size : max_sz;

    if (verbose) {
      std::print("threads requested: {}\n"
                 "reader threads: {}\n",
                 n_threads.workers, n_threads.readers);
      if (is_bgzf(input_format))
        std::print("decompression threads: {}\n", n_threads.decomp);
      std::println("memory requested: {}\n"
                   "tile analysis requested: {}\n"
                   "k-mer analysis requested: {}\n"
                   "use base groups in output: {}\n"
                   "input file format: {}\n"
                   "input files:",  //
                   size_to_units(buffer_size), do_tiles, do_kmers, do_groups,
                   format_description);
      std::ranges::for_each(infiles,
                            [](const auto &fn) { std::println("{}", fn); });
    }

    run_mode mode;
    mode.tiles(do_tiles);
    mode.kmers(do_kmers);
    mode.groups(do_groups);

    start_analysis(mode, buffer_size, n_threads, infos, infiles, outdirs);

    if (verbose)
      std::println(
        "total run time: {:.6g}s",
        duration(start_time, std::chrono::high_resolution_clock::now()));
  }
  catch (const std::exception &e) {
    std::println("{}", e.what());
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
