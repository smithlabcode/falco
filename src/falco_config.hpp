// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_FALCO_CONFIG_HPP_
#define SRC_FALCO_CONFIG_HPP_

#include <string>

class run_mode;

// ADS: run mode is an out-param because it might be partially set already
auto
load_config_and_set_graders(const std::string &filename,
                            run_mode &mode) -> void;

#endif  // SRC_FALCO_CONFIG_HPP_
