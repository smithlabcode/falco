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
    path_to_binary = std::string{path_buf};
#elif defined(_WIN32)
  const DWORD size =
    GetModuleFileName(nullptr, std::data(path_buf), path_buf_len);
  if (size > 0)
    path_to_binary = std::string{path_buf};
#else
  (void)path_buf;
#endif
  if (path_to_binary.empty())
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
