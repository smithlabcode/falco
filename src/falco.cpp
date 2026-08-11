// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

// clang-format off
static constexpr auto about = R"(Falco v{})";
static constexpr auto description =
  R"(Note: Always use bgzip when compressing files to gz format. It comes with samtools.

EXAMPLES:

Use all defaults making output files in 'results/SRX081761_1':
$ falco -o results SRX081761_1.fastq

Use 8 cores and analyze all files in 'project' with the fq suffix, and results
for each input file in their own subdirectory of 'results':
$ falco -t 8 -o results project/*.fq

Use 64 cores on multiple BAM files in the current directory:
$ falco -t 64 -o results sample*.bam

Skip adapter content and tile-specific analysis:
$ falco --no-adap --no-tiles -o results SRX081761.bam

Use configuration settings from 'my_limits.txt':
$ falco --config my_limits.txt -o results project_sample2.fq.gz

Default configuration files can be found here:
{}
Use these as templates. Copy and modify them to customize your analysis.
)";
// clang-format on

#include "adapter_set.hpp"
#include "bam_file.hpp"
#include "contaminants.hpp"
#include "falco_analyzer.hpp"
#include "falco_config.hpp"
#include "falco_file_format.hpp"
#include "falco_grade.hpp"
#include "falco_utils.hpp"
#include "fastq_bgzf_file.hpp"
#include "fastq_file.hpp"
#include "fastq_gz_file.hpp"
#include "get_binary_dir.hpp"
#include "original_duplicates.hpp"
#include "quality_score.hpp"
#include "reads_file.hpp"  // IWYU pragma: keep
#include "results_collector.hpp"
#include "results_summary.hpp"
#include "sam_file.hpp"
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

static auto
write_file(const auto &filename, const auto &data) {
  std::ofstream out(filename);
  if (!out)
    throw std::runtime_error("failed to open file: " + std::string(filename));
  std::print(out, "{}", data);
}

static auto
write_output(const run_mode &mode, std::vector<file_info> &infos,
             const std::vector<std::string> &outdirs,
             std::vector<results_collector> &&results) {
  static constexpr auto report_filename = "fastqc_data.txt";
  static constexpr auto html_filename = "fastqc_report.html";
  static constexpr auto summary_filename = "summary.txt";
  for (const auto [result, info, outdir] :
       std::views::zip(results, infos, outdirs)) {
    const auto outdir_path = std::filesystem::path{outdir};
    const auto summary = results_summary(std::move(result), mode, info);
    write_file(outdir_path / report_filename, summary.get_report());
    write_file(outdir_path / html_filename, summary.get_html());
    write_file(outdir_path / summary_filename, summary.get_summary());
  }
}

[[nodiscard]] static auto
make_reads_files(const std::vector<file_info> &infos,
                 const std::vector<std::string> &infiles,
                 const std::int64_t buf_size) -> std::vector<reads_file_t> {
  const auto n_infiles = std::size(infiles);
  std::vector<reads_file_t> reads_files;
  reads_files.reserve(n_infiles);
  for (const auto [infile, info] : std::views::zip(infiles, infos)) {
    if (info.format == falco::file_format::bam)
      reads_files.emplace_back(bam_file(infile, buf_size));
    else if (info.format == falco::file_format::sam)
      reads_files.emplace_back(sam_file(infile, buf_size));
    else if (info.format == falco::file_format::fastq_bgzf)
      reads_files.emplace_back(fastq_bgzf_file(infile, buf_size));
    else if (info.format == falco::file_format::fastq_gz)
      reads_files.emplace_back(fastq_gz_file(infile, buf_size));
    else if (info.format == falco::file_format::fastq)
      reads_files.emplace_back(fastq_file(infile, buf_size));
    else
      throw std::runtime_error("unsupported file format: " + info.description);
  }
  return reads_files;
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
      if (input_format == falco::file_format::sam)
        return estimate_n_reads_bam(infile);
      if (input_format == falco::file_format::fastq_bgzf)
        return estimate_n_reads_fastq_bgzf(infile);
      if (input_format == falco::file_format::fastq_gz)
        return estimate_n_reads_fastq_gz(infile);
      if (input_format == falco::file_format::fastq)
        return estimate_n_reads_fastq(infile);
      throw std::runtime_error("invalid reads file format");
    }();
    infos.push_back({
      .name = std::filesystem::path{infile}.filename().string(),
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

// ADS: in functions below, the input_record_multiplier is the relative space
// for a record: ideally for FASTQ it would be 2x the (max) size of the read
// name, and 2x the (max) sequence length. For BAM it would be 0.5x the sequence
// length, the name length and the length of the quality score sequence, plus
// all other BAM fields. The buffer_multiplier is how much total buffer space is
// needed for a given record, which is 3 for BAM and BGZF because one for
// decompression, and then two inflated buffers that are swapped so once can be
// analyzed while the other is the target of other inflation work.

[[nodiscard]] static auto
get_max_read_length(const auto buffer_size) {
  static constexpr auto input_record_multiplier = 2;
  static constexpr auto buffer_multiplier = 3;  // buffer capacity split
  static constexpr auto max_read_name_len = 64 * 1024;
  const auto max_read_len = (buffer_size - max_read_name_len) /
                            (buffer_multiplier * input_record_multiplier);
  return max_read_len;
}

[[nodiscard]] static auto
get_min_buffer_size(const auto max_read_length) {
  static constexpr auto input_record_multiplier = 2;
  static constexpr auto buffer_multiplier = 3;
  static constexpr auto max_read_name_len = 64 * 1024;
  const auto min_buffer_size =
    max_read_length * input_record_multiplier * buffer_multiplier +
    max_read_name_len;
  return min_buffer_size;
}

int
main(int argc, char *argv[]) {
  try {
    static constexpr auto buffer_size_default = 256 * 1024 * 1024;
    static constexpr std::int64_t min_buf_size = 1024 * 1024;
    std::vector<std::string> infiles;
    std::string contam_file;
    std::string config_file;
    std::string adapters_file;
    std::string outdir;
    std::int64_t buffer_size{buffer_size_default};

    int do_tiles{};
    int do_kmers{};
    int do_dup_analysis{};
    int do_adap{};

    int do_groups{};
    int do_bisulfite{};
    int do_original_dups{};

    std::uint32_t n_threads{1};
    std::uint32_t max_read_length{};

    int verbose{};

    using std::literals::string_literals::operator""s;
    const auto size_from_units = CLI::AsNumberWithUnit(std::map{
      std::pair{"G"s, gigabytes},
      {"M"s, megabytes},
      {"k"s, kilobytes},
    });

    const auto license_callback = [&](auto) {
      std::print("{}", license_text);
      throw CLI::Success();
    };

    struct FormatWithoutFlagDefaults : public CLI::Formatter {
      FormatWithoutFlagDefaults() : Formatter() {
        CLI::FormatterBase::enable_default_flag_values_ = false;
      }
    };

    CLI::App app{std::format(about, VERSION)};
    const auto fmt = std::make_shared<FormatWithoutFlagDefaults>();
    app.formatter(fmt);

    argv = app.ensure_utf8(argv);
    app.usage(
      std::format("Usage: {} [options] -o OUTDIR INFILES", PROJECT_NAME));
    if (argc >= 2)
      app.footer(std::format(description, falco::get_share_dir()));

    // NOLINTNEXTLINE (cppcoreguidelines-avoid-magic-numbers)
    app.get_formatter()->long_option_alignment_ratio(0.2);
    app.set_help_flag("-h,--help", "Print more detailed help");
    app.set_version_flag("--version", VERSION, "Print program version");
    // clang-format off
    app.add_flag("--license", license_callback, "Print full license")
      ->callback_priority(CLI::CallbackPriority::PreRequirementsCheck);
    app.add_option("INFILES", infiles,
                   "FASTQ (plain, GZIP or BGZF) or BAM/SAM")
      ->required()
      ->option_text(" ")
      ->check(CLI::ExistingFile);
    app.add_option("-o,--output", outdir, "Output directory")
      ->required()
      ->option_text("DIR");
    app.add_option("-t,--threads", n_threads,
                   std::format("Threads for analysis (this machine supports: {})",
                               std::thread::hardware_concurrency()))
      ->option_text(std::format("[{}]", n_threads));
    app.add_option("-m,--mem", buffer_size,
                   "Input memory buffer size (G/M/K units ok)")
      ->check(CLI::Range(min_buf_size, std::numeric_limits<std::int64_t>::max()))
      ->option_text(std::format("[{}]", size_to_units(buffer_size_default)))
      ->capture_default_str()
      ->transform(size_from_units);
    app.add_flag("-v,--verbose", verbose, "Print more info about the run")
      ->option_text(" ");
    app.add_option("--config", config_file,
                   "Config file (priority to command line options)")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("--contaminants", contam_file,
                   "File of non-default contaminant sequences to use")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("--adapters", adapters_file,
                   "File of non-default adapters sequences to use")
      ->option_text("FILE")
      ->check(CLI::ExistingFile);
    app.add_option("--max-length", max_read_length,
                   "Use for reads longer than 1M bp (G/M/K units ok)")
      ->option_text(" ")
      ->capture_default_str()
      ->transform(size_from_units);
    app.add_flag("--bisulfite", do_bisulfite,
                 "Assume bisulfite when grading sequence content")
      ->option_text(" ");
    const auto orig_dups_opt =
      app.add_flag("--orig-dups", [&](const auto x) {
        do_original_dups = x;
        do_dup_analysis = 1;
      }, "Use original duplication mode (turns dups on)")
      ->option_text(" ");
    app.add_flag("--groups", do_groups, "Group base positions in output")
      ->option_text(" ");
    app.add_flag("--tiles,!--no-tiles", do_tiles,
                 "Toggle per-tile quality analysis (default: on)")
      ->option_text(" ");
    app.add_flag("--dups,!--no-dups", do_dup_analysis,
                 "Toggle duplication/overrep analysis (default: on)")
      ->excludes(orig_dups_opt)
      ->option_text(" ");
    app.add_flag("--adap,!--no-adap", do_adap,
                 "Toggle adapter analysis (default: on)")
      ->option_text(" ");
    app.add_flag("--kmers,!--no-kmers", do_kmers,
                 "Toggle k-mer analysis (default: off)")
      ->option_text(" ");
    // clang-format on

    const auto start_time{std::chrono::high_resolution_clock::now()};

    if (argc < 2) {
      std::println("{}", app.help());
      return EXIT_SUCCESS;
    }
    CLI11_PARSE(app, argc, argv);

    run_mode mode;  // declare mode here so we can assign from config file
    if (!config_file.empty())
      load_config_and_set_graders(config_file, mode);
    else
      // if no config file, use default graders
      grader_set::instance();

    // now set run mode values to take priority over config file
    if (!adapters_file.empty())
      do_adap = 1;
    mode.set_do_adap(do_adap);
    mode.set_do_dups(do_dup_analysis);
    mode.set_do_overrep(do_dup_analysis);
    mode.set_do_kmers(do_kmers);
    mode.set_do_tiles(do_tiles);
    mode.set_do_groups(do_groups);
    mode.set_do_bisulfite(do_bisulfite);
    mode.set_do_original_dups(do_original_dups);
    mode.set_unassigned();

    const auto outdirs = make_outdirs(infiles, outdir);

    if (!contam_file.empty()) {
      load_contaminants(contam_file);
      if (verbose)
        std::print("contaminants file: {}\n"
                   "number of contaminants: {}\n",
                   contam_file, std::size(contaminants));
    }

    const auto &as = adapter_set::instance(mode, adapters_file);
    if (const auto [is_valid, message] = as.validate(); !is_valid) {
      std::println("{}", message);
      return EXIT_FAILURE;
    }

    if (!adapters_file.empty() && verbose)
      std::print("adapters file: {}\n"
                 "number of adapters: {}\n",
                 adapters_file, adapter_set::n_adapters());

    // not const because infos will change later when we can deduce the encoding
    auto infos = get_file_info(infiles);

    // restrict buffer size to avoid using a possibly harmful amount of memory
    const auto get_sz = [](const auto &i) { return i.size; };
    const auto max_sz = std::ranges::max(std::views::transform(infos, get_sz));
    buffer_size = buffer_size < max_sz ? buffer_size : max_sz;

    const auto min_buffer_size = get_min_buffer_size(max_read_length);
    if (min_buffer_size > buffer_size) {
      buffer_size = min_buffer_size;
      if (verbose)
        std::println("buffer size increased to accommodate max read length");
    }

    if (max_read_length == 0)
      max_read_length = get_max_read_length(buffer_size);

    if (verbose) {
      std::println("threads requested: {}\n"
                   "input memory buffer size: {}\n"
                   "max analyzable read length: {}\n"
                   "tile analysis requested: {}\n"
                   "k-mer analysis requested: {}\n"
                   "dups analysis requested: {}\n"
                   "adapter analysis requested: {}\n"
                   "use base groups in output: {}\n"
                   "input files:",
                   n_threads, size_to_units(buffer_size),
                   size_to_units(max_read_length, "bp"), mode.do_tiles(),
                   mode.do_kmers(), mode.do_dups(), mode.do_adap(),
                   mode.do_groups());
      std::ranges::for_each(infos, [](const auto &info) {
        std::println("{}\t{}\t{}", info.name, info.description,
                     size_to_units(info.size, std::string{}));
      });
    }

    std::vector<dups_init_t> dups =
      do_original_dups
        ? initialize_original_duplicates(infiles, infos, n_threads)
        : std::vector<dups_init_t>{};
    auto results =
      analyze(n_threads, mode, infos,
              make_reads_files(infos, infiles, buffer_size), std::move(dups));
    write_output(mode, infos, outdirs, std::move(results));

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
