// SPDX-License-Identifier: MIT; Copyright 2026 Andrew D Smith

#ifndef SRC_TASK_QUEUE_HPP_
#define SRC_TASK_QUEUE_HPP_

#include "falco_task.hpp"

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <iterator>
#include <list>
#include <mutex>
#include <utility>
#include <variant>

struct task_queue {
  std::mutex mtx;
  std::condition_variable cv;
  std::list<std::pair<std::int32_t, task_t>> q;  // file_id and task
  bool stop{};

  auto
  push(const auto file_id, task_t task) {
    std::unique_lock tq_lock(mtx);
    q.emplace_back(file_id, std::move(task));
    tq_lock.unlock();
    cv.notify_one();
  }

  [[nodiscard]] auto
  wait_and_acquire_lock() -> std::unique_lock<std::mutex> {
    std::unique_lock lock(mtx);
    cv.wait(lock, [this] { return !q.empty() || stop; });
    return lock;
  }

  [[nodiscard]] auto
  pop_and_release_lock(auto &tq_lock) {
    auto task = std::move(q.front());
    q.pop_front();
    tq_lock.unlock();
    cv.notify_one();
    return task;
  }

  [[nodiscard]] auto
  is_finished() const -> bool {
    return stop && q.empty();
  }

  auto
  request_shutdown() {
    stop = true;
    cv.notify_all();
  }
};

#endif  // SRC_TASK_QUEUE_HPP_
