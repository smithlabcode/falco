// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "sam_file.hpp"
#include "bamrec.hpp"      // IWYU pragma: keep
#include "task_queue.hpp"  // IWYU pragma: keep

#include <htslib/hfile.h>  // IWYU pragma: keep
#include <htslib/sam.h>    // IWYU pragma: keep

#include <atomic>
#include <cstdint>

auto
sam_file::load_next([[maybe_unused]] const std::int32_t file_id,   //
                    [[maybe_unused]] task_queue &tq,               //
                    [[maybe_unused]] std::atomic_int32_t &n_tasks  //
                    ) -> void {}

auto
sam_file::get_chunks([[maybe_unused]] const std::int64_t n_chunks,  //
                     [[maybe_unused]] const std::int32_t file_id,   //
                     [[maybe_unused]] task_queue &tq,               //
                     [[maybe_unused]] std::atomic_int32_t &n_tasks  //
                     ) -> void {}

auto
sam_file::make_tasks(const std::int64_t n_chunks, const std::int32_t file_id,
                     task_queue &tq, std::atomic_int32_t &n_tasks) -> void {
  get_chunks(n_chunks, file_id, tq, n_tasks);
}
