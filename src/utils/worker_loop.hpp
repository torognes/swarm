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

#ifndef SWARM_UTILS_WORKER_LOOP_H
#define SWARM_UTILS_WORKER_LOOP_H

#include <mutex>


/*
  run_worker_loop drives one worker thread's self-scheduling loop, the
  pattern shared by swarm's dynamically load-balanced ThreadRunner
  passes (ported from vsearch, whose ThreadRunner is this one's twin):
  repeatedly claim the next unit of work under a shared input mutex,
  then process that unit with the mutex released so the other workers
  can claim in parallel.

  claim is called while input_mutex is held. It should claim the next
  unit of work (advance the shared cursor, and record which unit was
  claimed somewhere the work callable can see it) and return true, or
  return false when the input is exhausted. Returning false ends the
  loop. A claim that has to skip over units no worker should process
  simply keeps advancing the cursor until it can claim or report
  exhaustion.

  work is called with no lock held. It performs the computation for the
  unit the preceding claim call selected, communicating through state
  the two callables share (typically captured by reference). Any output
  synchronisation is work's own responsibility -- taking input_mutex
  again for a final accumulation step is fine, since the loop holds no
  lock while work runs.

  Centralising the loop makes "release the input mutex before doing the
  work" a property of one place instead of an invariant each caller
  re-establishes by hand-placing input_lock.unlock().

  This cannot be marked noexcept: locking the mutex may throw
  std::system_error and the supplied callables are unconstrained. This
  matches ThreadRunner, whose run() is likewise not noexcept.
*/
template <typename ClaimFn, typename WorkFn>
auto run_worker_loop(std::mutex & input_mutex,
                     ClaimFn claim,
                     WorkFn work) -> void {
  while (true) {
    {
      std::lock_guard<std::mutex> const lock(input_mutex);
      if (not claim()) {
        break;
      }
    }
    work();
  }
}

#endif  // SWARM_UTILS_WORKER_LOOP_H
