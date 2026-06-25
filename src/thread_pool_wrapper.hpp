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

#ifndef SRC_THREAD_POOL_WRAPPER_HPP_
#define SRC_THREAD_POOL_WRAPPER_HPP_

#include <htslib/sam.h>
#include <htslib/thread_pool.h>

#include <cstdint>
#include <stdexcept>

[[nodiscard]] static auto
make_thread_pool(const auto n_threads) {
  return n_threads > 0 ? hts_tpool_init(n_threads) : nullptr;
}

struct falco_thread_pool {
  htsThreadPool t{};

  explicit falco_thread_pool(const std::uint64_t n_threads) :
    t{make_thread_pool(n_threads), 0} {
    if (n_threads > 0 && t.pool == nullptr)
      throw std::runtime_error("failed to construct thread pool");
  }

  [[nodiscard]] auto
  n_threads() const -> std::uint64_t {
    return t.pool ? hts_tpool_size(t.pool) : 0LU;
  }

  // clang-format off
  falco_thread_pool(const falco_thread_pool &) = delete;
  auto operator=(const falco_thread_pool &) -> falco_thread_pool & = delete;
  falco_thread_pool(falco_thread_pool &&) noexcept = delete;
  auto operator=(falco_thread_pool &&) noexcept -> falco_thread_pool & = delete;
  ~falco_thread_pool() { if (t.pool) hts_tpool_destroy(t.pool); }
  // clang-format on
};

#endif  // SRC_THREAD_POOL_WRAPPER_HPP_
