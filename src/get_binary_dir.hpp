// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_GET_BINARY_DIR_HPP_
#define SRC_GET_BINARY_DIR_HPP_

#include <string>

namespace falco {

[[nodiscard]] auto
get_share_dir() -> std::string;

[[nodiscard]] auto
get_binary_dir() -> std::string;

}  // namespace falco

#endif  // SRC_GET_BINARY_DIR_HPP_
