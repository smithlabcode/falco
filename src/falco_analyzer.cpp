// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "falco_analyzer.hpp"

#include "bamrec.hpp"
#include "bgzf_block.hpp"
#include "duplication_results.hpp"
#include "falco_utils.hpp"
#include "file_info.hpp"
#include "fqrec.hpp"
#include "reads_file.hpp"
#include "results_collector.hpp"
#include "samrec.hpp"
#include "task_queue.hpp"

#include <atomic>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <mutex>
#include <ranges>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

class run_mode;

[[nodiscard]] auto
analyze(const std::uint32_t n_threads, const run_mode &mode,
        std::vector<file_info> &infos, std::vector<reads_file_t> reads_files,
        std::vector<dups_init_t> dups_init) -> std::vector<results_collector> {
  assert(std::size(reads_files) == std::size(infos));
  const std::int32_t n_files = static_cast<std::int32_t>(std::size(infos));
  if (dups_init.empty())
    dups_init.resize(n_files);
  std::vector<std::atomic_int32_t> n_tasks(n_files);
  std::atomic_uint32_t n_active_files{static_cast<std::uint32_t>(n_files)};
  auto results =
    std::vector(n_threads, std::vector<results_collector>(n_files));

  // add initial jobs so workers can work immediately
  task_queue tq;
  for (const auto file_id : std::views::iota(0, n_files))
    tq.push(file_id, std::monostate{});

  {  // scope to join jthreads
    std::vector<std::jthread> workers;
    for (const auto th_id : std::views::iota(0u, n_threads))
      workers.emplace_back([&, n_threads, th_id] {
        auto &res = results[th_id];
        for (const auto [file_id, info] : falco::views::enumerate(infos))
          res[file_id].init(mode, info, dups_init[file_id]);
        while (true) {
          auto tq_lock = tq.wait_and_acquire_lock();
          if (tq.is_finished())
            return;
          auto [file_id, task] = tq.pop_and_release_lock(tq_lock);
          if (std::holds_alternative<fq_task_t>(task))
            process_reads(res[file_id], std::get<fq_task_t>(task));
          else if (std::holds_alternative<bam_task_t>(task))
            process_reads(res[file_id], std::get<bam_task_t>(task));
          else if (std::holds_alternative<sam_task_t>(task))
            process_reads(res[file_id], std::get<sam_task_t>(task));
          else if (std::holds_alternative<bgzf_block_t>(task))
            decompress(std::get<bgzf_block_t>(task));
          else {  // monostate means read more data
            if (!is_active(reads_files[file_id])) {
              assert(n_tasks[file_id] == 0);
              reset(reads_files[file_id]);
              if (n_active_files.fetch_sub(1, std::memory_order_relaxed) == 1) {
                // ADS: it would not be wrong it two different threads arrived
                // here, but using fetch_sub should prevent that anyway.
                tq.request_shutdown();
                return;
              }
              // Must 'continue' here to avoid decrementing n_tasks[file_id]
              // below, since within this branch we might exit before
              // decrementing to remove the value that we would have added; the
              // same is not true of the 'else' case below when make_tasks is
              // called. So the ++n_tasks[file_id] is inside that function in
              // each case of typeof(reads_files). Unfortunately changes to
              // n_tasks need to happen in different places.
              continue;
            }
            else
              make_tasks(reads_files[file_id], n_threads, file_id, tq,
                         n_tasks[file_id]);
          }
          if (n_tasks[file_id].fetch_sub(1, std::memory_order_relaxed) == 1) {
            assert(n_tasks[file_id] == 0);
            tq.push(file_id, std::monostate{});
          }
        }
      });
  }

  std::vector<results_collector> finalized;
  finalized.reserve(n_files);
  for (const auto file_id : std::views::iota(0, n_files)) {
    accumulate_results(results, file_id);
    finalized.emplace_back(std::move(results.front()[file_id]));
    finalized.back().finalize(infos[file_id]);
  }

  return finalized;
}
