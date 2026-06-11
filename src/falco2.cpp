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

#include "bam_file.hpp"
#include "contaminants.hpp"
#include "falco_analyzer.hpp"
#include "falco_file_format.hpp"
#include "falco_results.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "quality_score.hpp"
#include "tile_processor.hpp"

#include "format_output.hpp"

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

template <typename results_t>
static auto
run(auto &infos, auto &reads_files, const auto n_threads,
    const auto &outfiles) {
  using rec_t = std::decay_t<decltype(reads_files)>::value_type::rec_t;
  const auto n_files = std::size(reads_files);

  analyzer_t<results_t, rec_t> analyzer(n_threads, n_files, reads_files, infos);

  for (const auto file_id : std::views::iota(0u, n_files))
    accumulate_results(analyzer.results, file_id);

  for (const auto [results, info, outfile] :
       std::views::zip(analyzer.results.front(), infos, outfiles)) {
    set_quality_score_encoding(results.qual_by_pos, info);
    if (!is_mapped_reads(info.format))
      results.finalize_qual_encoding(info.encoding);
    grades g;
    const auto report = results.get_report(info, g);

    std::ofstream out(outfile);
    if (!out)
      throw std::runtime_error("failed to open file: " + outfile);
    std::print(out, "{}", report);
  }
}

static auto
run_mode_selector(const run_mode mode, std::vector<file_info> &infos,
                  auto &reads_files, const auto n_threads,
                  const auto &outfiles) {
  if (tiles(mode) && kmers(mode))
    run<falco_results_tile_kmer>(infos, reads_files, n_threads, outfiles);
  else if (tiles(mode))
    run<falco_results_tile>(infos, reads_files, n_threads, outfiles);
  else if (kmers(mode))
    run<falco_results_kmer>(infos, reads_files, n_threads, outfiles);
  else
    run<falco_results>(infos, reads_files, n_threads, outfiles);
}

static auto
start_analysis(const run_mode mode, const auto buf_size, const auto n_threads,
               const auto input_format, const auto &format_description,
               auto &infos, const auto &infiles, const auto &outfiles) {
  // ADS: this function is bloated mostly to avoid allowing default construction
  // or move assignment in the fastq_file, fastq_gz_file, etc.
  const auto n_infiles = std::size(infiles);
  if (is_mapped_reads(input_format)) {
    falco_thread_pool p(n_threads > 1 ? n_threads - 1 : 1);
    std::vector<bam_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size, p); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outfiles);
  }
  else if (input_format == falco::file_format::fastq_bgzf) {
    falco_thread_pool p(n_threads > 1 ? n_threads - 1 : 1);
    std::vector<fastq_bgzf_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size, p); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outfiles);
  }
  else if (input_format == falco::file_format::fastq_gz) {
    std::vector<fastq_gz_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outfiles);
  }
  else if (input_format == falco::file_format::fastq) {
    std::vector<fastq_file> f;
    f.reserve(n_infiles);
    const auto m = [&](const auto &i) { f.emplace_back(i, buf_size); };
    std::ranges::for_each(infiles, m);
    run_mode_selector(mode, infos, f, n_threads, outfiles);
  }
  else {
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
    const auto [n_reads_est, filesize] = [&] {
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
    static constexpr auto outfile_name = "fastqc_data.txt";
    static constexpr auto buf_size_default = 512 * 1024 * 1024;
    std::vector<std::string> infiles;
    std::string contam_file;
    std::string outdir;
    std::int64_t buf_size{buf_size_default};
    std::uint32_t n_threads{1};
    std::uint32_t n_bam_threads{1};

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
                   "Input file: FASTQ (plain, gz or bgzf) or BAM/SAM")
      ->required()
      ->option_text("FILES")
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outdir, "Output directory")
      ->required()
      ->option_text("FILE");
    app.add_option("-c,--contaminants", contam_file,
                   "File of contaminant sequences to use")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("-b,--bam-threads", n_bam_threads,
                   "Threads to use for decompression")
      ->option_text(std::format("[{}]", n_bam_threads));
    app.add_option("-t,--threads", n_threads,
                   std::format("Threads to use (this machine supports: {})",
                               std::thread::hardware_concurrency()))
      ->option_text(std::format("[{}]", n_threads));
    app.add_option("-m,--mem", buf_size,
                   "Memory buffer size for IO (G/M/K units ok)")
      ->option_text(std::format("[{}]", size_to_units(buf_size_default)))
      ->capture_default_str()
      ->transform(size_from_units);
    app.add_flag("-v,--verbose", verbose, "Print more info while running.");
    app.add_flag("--tiles", do_tiles, "Enable per-tile analysis");
    app.add_flag("--kmers", do_kmers, "Enable k-mer analysis");
    // clang-format on

    const auto start_time{std::chrono::high_resolution_clock::now()};

    if (argc < 2) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    const auto outdirs = make_outdirs(infiles, outdir);

    if (!contam_file.empty()) {
      load_contaminants(contam_file);
      if (verbose)
        std::print("contaminants file: {}\n"
                   "number of contaminants: {}\n",
                   contam_file, std::size(contaminants));
    }

    if (n_threads > std::thread::hardware_concurrency())
      n_threads = std::thread::hardware_concurrency();

    if (verbose)
      std::print("threads requested: {}\n"
                 "memory requested: {}\n"
                 "tile analysis requested: {}\n"
                 "k-mer analysis requested: {}\n",
                 n_threads, size_to_units(buf_size), do_tiles, do_kmers);

    auto infos = get_file_info(infiles);

    const auto format_description = infos.front().description;
    const auto input_format = infos.front().format;
    const auto ok_fmt = [&](const auto &x) { return x.format == input_format; };
    if (!std::ranges::all_of(infos, ok_fmt)) {
      std::println("unequal input formats [expected {}]", format_description);
      return EXIT_FAILURE;
    }

    // restrict buffer size to avoid using a possibly harmful amount of memory
    const auto get_sz = [](const auto &i) { return i.size; };
    const auto max_sz = std::ranges::max(std::views::transform(infos, get_sz));
    buf_size = buf_size < max_sz ? buf_size : max_sz;

    if (verbose) {
      std::println("input file format: {}", format_description);
      std::println("input files:");
      for (auto infile : infiles)
        std::println("{}", infile);
    }

    run_mode mode;
    mode.tiles(do_tiles);
    mode.kmers(do_kmers);

    const auto outfiles = std::views::transform(outdirs, [](const auto &d) {
      return (std::filesystem::path{d} / outfile_name).string();
    });

    start_analysis(mode, buf_size, n_threads, input_format, format_description,
                   infos, infiles, outfiles);

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
