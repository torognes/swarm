/*
    SWARM

    Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>


class ThreadRunner {
private:

  struct thread_s {
    int64_t thread_id {0};
    std::function<void(int64_t)> fun;
    std::thread thread;
    std::mutex workmutex;
    std::condition_variable workcond;
    int64_t work {0}; /* 1: work available, 0: wait, -1: quit */
  };

  std::vector<struct thread_s> thread_array;

  static auto worker(struct thread_s * tip) -> void {
    std::unique_lock<std::mutex> lock(tip->workmutex);

    /* loop until signalled to quit */
    while (tip->work >= 0) {
      /* wait for work available */
      if (tip->work == 0) {
        tip->workcond.wait(lock);
      }

      if (tip->work > 0) {
        tip->fun(tip->thread_id);
        tip->work = 0;
        tip->workcond.notify_one();
      }
    }
  }


public:

  // refactoring: heaptrack detects a memory leak of 640 bytes for
  // each thread created by this member function. Backtrace:
  // ThreadRunner::ThreadRunner(int, std::function<void(long)> const&)
  //   __pthread_create_2_1 in libc.so.6
  //   allocate_stack in libc.so.6
  //   __GI__dl_allocate_tls in ld-linux-x86-64.so.2
  //   allocate_dtv in ld-linux-x86-64.so.2
  //   calloc in ld-linux-x86-64.so.2
  ThreadRunner(int thread_count,
               const std::function<void(int64_t nth_thread)> & function) :
      thread_array(static_cast<std::size_t>(thread_count)) {
    /* init and create worker threads */
    auto counter = 0LL;
    for (auto & tip: thread_array) {
        tip.thread_id = counter;
        tip.fun = function;
        tip.thread = std::thread(worker, &tip);
        ++counter;
      }
  }


  ~ThreadRunner() {
    /* ask threads to quit */
    /* wait for them to join */

    for (auto & tip: thread_array) {
        /* tell worker to quit */
        {
          const std::lock_guard<std::mutex> lock(tip.workmutex);
          tip.work = -1;
          tip.workcond.notify_one();
        }
        /* wait for worker to quit */
        tip.thread.join();
    }
  }

  ThreadRunner(const ThreadRunner&) = delete; // copy constructor
  ThreadRunner(ThreadRunner&&) = delete; // move constructor
  auto operator=(const ThreadRunner&) -> ThreadRunner& = delete; // copy assignment constructor
  auto operator=(ThreadRunner&&) -> ThreadRunner& = delete; // move assignment constructor

  auto run() -> void {
    /* wake up threads */
    for (auto & tip: thread_array) {
        const std::lock_guard<std::mutex> lock(tip.workmutex);
        tip.work = 1;
        tip.workcond.notify_one();
    }

    /* wait for threads to finish their work */
    for (auto & tip: thread_array) {
        std::unique_lock<std::mutex> lock(tip.workmutex);
        tip.workcond.wait(lock, [&tip]() -> bool { return tip.work <= 0; });
    }
  }
};
