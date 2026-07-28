// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_TASK_HPP_
#define SRC_FALCO_TASK_HPP_

#include "bamrec.hpp"
#include "bgzf_block.hpp"
#include "fqrec.hpp"
#include "samrec.hpp"

#include <variant>

using task_t =
  std::variant<std::monostate, fq_task_t, bam_task_t, sam_task_t, bgzf_block_t>;

#endif  // SRC_FALCO_TASK_HPP_
