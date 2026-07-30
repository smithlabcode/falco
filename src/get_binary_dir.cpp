// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#include "get_binary_dir.hpp"

#include <config.h>

#if defined(__linux__)
#include <limits.h>  // IWYU pragma: keep
#include <unistd.h>
#elif defined(__APPLE__)
#include <libproc.h>
#include <unistd.h>
#elif defined(_WIN32)
#include <iostream>
#include <windows.h>
#endif

#include <sys/types.h>  // for ssize_t

#include <array>
#include <filesystem>
#include <iterator>
#include <string>

namespace falco {

[[nodiscard]] auto
get_binary_dir() -> std::string {
  static constexpr auto path_buf_len = 1024;
  std::array<char, path_buf_len> path_buf{};

  std::string path_to_binary;
#if defined(__linux__)
  static constexpr auto exe_path = "/proc/self/exe";
  const ssize_t length =
    readlink(exe_path, std::data(path_buf), path_buf_len - 1);
  if (length != -1)
    path_to_binary = std::data(path_buf);
#elif defined(__APPLE__)
  const pid_t pid = getpid();
  const ssize_t length = proc_pidpath(pid, std::data(path_buf), path_buf_len);
  if (length > 0)
    path_to_binary = std::data(path_buf);
#elif defined(_WIN32)
  const DWORD size =
    GetModuleFileName(nullptr, std::data(path_buf), path_buf_len);
  if (size > 0)
    path_to_binary = std::data(path_buf);
#else
  (void)path_buf;
#endif
  if (path_to_binary.empty())  // cppcheck-suppress knownConditionTrueFalse
    return std::string{};
  return std::filesystem::path{path_to_binary}.parent_path();
}

[[nodiscard]] auto
get_share_dir() -> std::string {
  static constexpr auto bin_name = "bin";
  auto binary_dir = std::filesystem::path{get_binary_dir()};
  if (binary_dir.filename() == bin_name)
    binary_dir = binary_dir.parent_path();
  const auto share_path = binary_dir / DATADIR / PROJECT_NAME;
  return is_directory(share_path) ? share_path.string() : std::string{};
}

}  // namespace falco
