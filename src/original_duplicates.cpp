// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "original_duplicates.hpp"

#include "bam_file.hpp"
#include "duplication_results.hpp"
#include "falco_file_format.hpp"
#include "fastq_bgzf_file.hpp"
#include "file_info.hpp"

#include "boost/boost_unordered.hpp"

#include <algorithm>
#include <cassert>
#include <iterator>
#include <mutex>
#include <ranges>
#include <string>
#include <thread>
#include <vector>

[[nodiscard]] auto
initialize_original_duplicates(
  const std::vector<std::string> &infiles,
  [[maybe_unused]] const std::vector<file_info> &infos,
  [[maybe_unused]] const std::uint32_t n_threads) -> std::vector<dups_init_t> {
  assert(std::size(infiles) == std::size(infos));
  const auto n_files = std::size(infiles);
  std::vector<dups_map_t> dups(n_files);
  auto n_active_files = n_files;
  std::vector<std::jthread> workers;
  const auto n_workers =
    std::min(n_files, static_cast<std::uint64_t>(n_threads));
  workers.reserve(n_workers);
  std::mutex mtx;
  // NOLINTNEXTLINE(clang-analyzer-deadcode.DeadStores)
  for (const auto _ : std::views::iota(0u, n_workers))
    workers.emplace_back([&] {
      std::uint64_t file_id{};
      while (true) {
        {
          std::scoped_lock l(mtx);
          if (n_active_files == 0)
            return;
          file_id = --n_active_files;
        }
        dups[file_id] =
          falco::is_mapped_reads(infos[file_id].format)
            ? init_dups(infiles[file_id],
                        duplication_results::max_reads_to_hash)
            : init_dups_fq(infiles[file_id],
                           duplication_results::max_reads_to_hash);
      }
    });
  std::ranges::for_each(workers, [](auto &w) { w.join(); });
  std::vector<dups_init_t> ret;
  ret.reserve(n_files);
  for (const auto &d : dups)
    ret.emplace_back(d);
  return ret;
}
