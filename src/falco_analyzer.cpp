// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "falco_analyzer.hpp"

#include "bam_file.hpp"
#include "falco_task.hpp"
#include "falco_utils.hpp"
#include "fastq_file.hpp"
#include "file_info.hpp"
#include "reads_file.hpp"
#include "results_collector.hpp"
#include "run_mode.hpp"
#include "task_queue.hpp"

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <mutex>
#include <ranges>
#include <string>
#include <thread>
#include <utility>
#include <variant>
#include <vector>

analyzer_t::analyzer_t(const std::uint32_t n_threads, const run_mode &mode,
                       const std::vector<file_info> &infos,
                       std::vector<reads_file_t> &reads_files) :
  n_tasks(static_cast<std::int32_t>(std::size(infos))),
  n_active_files{static_cast<std::uint32_t>(std::size(infos))},
  results(n_threads, std::vector<results_collector>(std::size(infos))) {
  assert(std::size(reads_files) == std::size(infos));

  const std::uint32_t n_files = std::size(infos);

  const auto n_chunks = n_threads * n_chunks_per_thread;
  // set per-file information used to do the analysis
  for (auto &res : results)
    for (const auto [file_id, info] : std::views::enumerate(infos))
      res[file_id].init(mode, info);

  // add initial jobs so workers can work immediately
  for (const auto file_id : std::views::iota(0u, n_files))
    tq.push(file_id, std::monostate{});

  std::vector<std::jthread> workers;
  workers.reserve(n_threads);
  for (const auto th_id : std::views::iota(0u, n_threads))
    workers.emplace_back([this, n_chunks, th_id, &reads_files] {
      auto &res = results[th_id];
      while (true) {
        auto tq_lock = tq.wait_and_acquire_lock();
        if (tq.is_finished())
          return;
        auto [file_id, task] = tq.pop_and_release_lock(tq_lock);
        if (std::holds_alternative<fq_task_t>(task))
          process_reads(res[file_id], std::move(std::get<fq_task_t>(task)));
        else if (std::holds_alternative<bam_task_t>(task))
          process_reads(res[file_id], std::move(std::get<bam_task_t>(task)));
        else if (std::holds_alternative<sam_task_t>(task))
          process_reads(res[file_id], std::move(std::get<sam_task_t>(task)));
        else if (std::holds_alternative<bgzf_block_t>(task))
          decompress(std::move(std::get<bgzf_block_t>(task)));
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
            // below, since within this branch we might exit before decrementing
            // to remove the value that we would have added; the same is not
            // true of the 'else' case below when make_tasks is called. So the
            // ++n_tasks[file_id] is inside that function in each case of
            // typeof(reads_files). Unfortunately changes to n_tasks need to
            // happen in different places.
            continue;
          }
          else
            make_tasks(reads_files[file_id], n_chunks, file_id, tq,
                       n_tasks[file_id]);
        }
        if (n_tasks[file_id].fetch_sub(1, std::memory_order_relaxed) == 1) {
          assert(n_tasks[file_id] == 0);
          tq.push(file_id, std::monostate{});
        }
      }
    });
}
