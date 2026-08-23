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

#ifndef SWARM_UTILS_THREADS_H
#define SWARM_UTILS_THREADS_H

#include <algorithm>  // std::for_each
#include <cassert>
#include <condition_variable>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>  // std::next
#include <mutex>
#include <thread>
#include <vector>


class ThreadRunner {
private:

  enum struct Work_state : std::int8_t { wait = 0, work = 1, quit = -1 };

  struct thread_s {
    uint64_t thread_id {0};
    std::function<void(uint64_t)> fun;
    std::thread thread;
    std::mutex workmutex;
    std::condition_variable workcond;
    Work_state work {Work_state::wait};
  };

  std::vector<struct thread_s> thread_array;

  static auto worker(struct thread_s & tip) -> void {
    std::unique_lock<std::mutex> lock(tip.workmutex);

    /* loop until signalled to quit */
    while (tip.work != Work_state::quit) {
      /* wait for work available */
      if (tip.work == Work_state::wait) {
        tip.workcond.wait(lock);
      }

      if (tip.work == Work_state::work) {
        tip.fun(tip.thread_id);
        tip.work = Work_state::wait;
        tip.workcond.notify_one();
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
  ThreadRunner(std::size_t const thread_count,
               std::function<void(uint64_t nth_thread)> const & function) :
      thread_array(thread_count) {
    /* init and create worker threads */
    // std::ref is required: std::thread decays its arguments by
    // default, so passing `tip` directly would copy thread_s — which
    // contains non-copyable members (std::mutex, std::condition_variable,
    // std::thread) and would fail to compile. std::ref preserves the
    // reference semantics so worker() receives the live thread_s.
    uint64_t counter {0};
    for (auto & tip: thread_array) {
        tip.thread_id = counter;
        tip.fun = function;
        tip.thread = std::thread(worker, std::ref(tip));
        ++counter;
      }
  }


  ~ThreadRunner() {
    /* ask threads to quit */
    /* wait for them to join */
    for (auto & tip: thread_array) {
        /* tell worker to quit */
        {
          std::lock_guard<std::mutex> const lock(tip.workmutex);
          tip.work = Work_state::quit;
          tip.workcond.notify_one();
        }
        /* wait for worker to quit */
        tip.thread.join();
    }
  }


  ThreadRunner(ThreadRunner const &) = delete; // copy constructor
  ThreadRunner(ThreadRunner&&) = delete; // move constructor
  auto operator=(ThreadRunner const &) -> ThreadRunner& = delete; // copy assignment constructor
  auto operator=(ThreadRunner&&) -> ThreadRunner& = delete; // move assignment constructor


  // Wake the first 'first_n' workers only, and wait for those alone. A
  // caller with less work than configured threads asks for the count it can
  // keep busy and pays the wake-and-join latency for that count, not for the
  // full pool. The workers left out keep sleeping and their state is not
  // touched, so they must not be counted on to have current work to do.
  auto run(std::size_t const first_n) -> void {
    assert(first_n >= 1);
    assert(first_n <= thread_array.size());
    auto const past_last = std::next(thread_array.begin(),
                                     static_cast<std::ptrdiff_t>(first_n));

    /* wake up threads */
    std::for_each(thread_array.begin(), past_last,
                  [](struct thread_s & tip) -> void {
                    std::lock_guard<std::mutex> const lock(tip.workmutex);
                    tip.work = Work_state::work;
                    tip.workcond.notify_one();
                  });

    /* wait for threads to finish their work */
    std::for_each(thread_array.begin(), past_last,
                  [](struct thread_s & tip) -> void {
                    std::unique_lock<std::mutex> lock(tip.workmutex);
                    tip.workcond.wait(lock, [&tip]() -> bool { return tip.work != Work_state::work; });
                  });
  }

  auto run() -> void {
    run(thread_array.size());
  }
};

#endif
