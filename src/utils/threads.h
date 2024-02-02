/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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

#include <cstdint>
#include <vector>
#include <pthread.h>  // refactoring: C++11 replace with std::thread
#include "fatal.h"


class ThreadRunner
{
private:

  pthread_attr_t attr {};

  struct thread_s
  {
    int64_t thread_id;
    void (*fun)(int64_t thread_id);
    pthread_t pthread;
    pthread_mutex_t workmutex;
    pthread_cond_t workcond;
    int64_t work; /* 1: work available, 0: wait, -1: quit */
  };

  std::vector<struct thread_s> thread_array;

  static auto worker(void * void_ptr) -> void *
  {
    auto * tip = static_cast<struct thread_s *>(void_ptr);

    pthread_mutex_lock(&tip->workmutex);  // refactoring: prefer mutex_lockguard (RAII)

    /* loop until signalled to quit */
    while (tip->work >= 0)
      {
        /* wait for work available */
        if (tip->work == 0) {
          pthread_cond_wait(&tip->workcond, &tip->workmutex);
        }

        if (tip->work > 0)
          {
            (*tip->fun)(tip->thread_id);
            tip->work = 0;
            pthread_cond_signal(&tip->workcond);
          }
      }

    pthread_mutex_unlock(&tip->workmutex);

    return nullptr;
  }


public:

  // refactoring: heaptrack detects a memory leak of 640 bytes for
  // each thread created by this member function. Backtrace:
  // ThreadRunner::ThreadRunner(int, void (*)(long))
  //   __pthread_create_2_1 in libc.so.6
  //   allocate_stack in libc.so.6
  //   __GI__dl_allocate_tls in ld-linux-x86-64.so.2
  //   allocate_dtv in ld-linux-x86-64.so.2
  //   calloc in ld-linux-x86-64.so.2
  ThreadRunner(int thread_count, void (*function_ptr)(int64_t nth_thread))
  {

    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    /* allocate memory for thread data */
    thread_array.resize(static_cast<uint64_t>(thread_count));

    /* init and create worker threads */
    auto counter = 0LL;
    for(auto& tip: thread_array) {
        tip.thread_id = counter;
        tip.work = 0;
        tip.fun = function_ptr;
        pthread_mutex_init(&tip.workmutex, nullptr);
        pthread_cond_init(&tip.workcond, nullptr);
        if (pthread_create(&tip.pthread,
                           &attr,
                           worker,
                           static_cast<void*>(&tip)) != 0) {
          fatal(error_prefix, "Cannot create thread.");
        }
        ++counter;
      }
  }


  ~ThreadRunner()
  {
    /* ask threads to quit */
    /* sleep until they have quit */
    /* destroy threads */
    /* finish and clean up worker threads */

    for(auto& tip: thread_array) {
        /* tell worker to quit */
        pthread_mutex_lock(&tip.workmutex);
        tip.work = -1;
        pthread_cond_signal(&tip.workcond);
        pthread_mutex_unlock(&tip.workmutex);

        /* wait for worker to quit */
        if (pthread_join(tip.pthread, nullptr) != 0) {
          fatal(error_prefix, "Cannot join thread.");
        }

        pthread_cond_destroy(&tip.workcond);
        pthread_mutex_destroy(&tip.workmutex);
    }

    pthread_attr_destroy(&attr);
  }

  ThreadRunner(const ThreadRunner&) = delete; // copy constructor
  ThreadRunner(ThreadRunner&&) = delete; // move constructor
  auto operator=(const ThreadRunner&) -> ThreadRunner& = delete; // copy assignment constructor
  auto operator=(ThreadRunner&&) -> ThreadRunner& = delete; // move assignment constructor

  auto run() -> void {
    /* wake up threads */
    for(auto& tip: thread_array) {
        pthread_mutex_lock(&tip.workmutex);
        tip.work = 1;
        pthread_cond_signal(&tip.workcond);
        pthread_mutex_unlock(&tip.workmutex);
    }

    /* wait for threads to finish their work */
    for(auto& tip: thread_array) {
        pthread_mutex_lock(&tip.workmutex);
        while (tip.work > 0) {
          pthread_cond_wait(&tip.workcond, &tip.workmutex);
        }
        pthread_mutex_unlock(&tip.workmutex);
    }
  }
};
