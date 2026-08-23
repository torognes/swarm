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

#ifndef SWARM_UTILS_QGRAM_H
#define SWARM_UTILS_QGRAM_H


#include "cpu_features.hpp"  // Cpu_features
#include "algo_internal.hpp"  // ampliconinfo_s
#include "qgram_array.hpp"  // Qgram_store
#include "span.hpp"  // Span<uint64_t>
#include "thread_count.hpp"  // ThreadCount
#include "threads.hpp"  // ThreadRunner
#include "view.hpp"  // View<uint64_t>
#include <cstdint>  // uint64_t
#include <vector>


struct Parameters;  // defined in swarm.hpp
class Data;          // defined in db.hpp


// RAII wrapper around the per-seed qgram-distance dispatch: builds the
// per-sequence qgram store, owns it together with the worker
// ThreadRunner, a per-thread scratch vector, and the cached
// Cpu_features. The destructor joins the worker threads — no explicit
// teardown call needed.
class QgramDiffer {
public:
  QgramDiffer(struct Parameters const & parameters,
              Data const & data);

  // Non-copyable, non-movable: the ThreadRunner's lambda captures
  // `this`, so the object must keep a stable address.
  QgramDiffer(QgramDiffer const &) = delete;
  QgramDiffer(QgramDiffer &&) = delete;
  auto operator=(QgramDiffer const &) -> QgramDiffer & = delete;
  auto operator=(QgramDiffer &&) -> QgramDiffer & = delete;
  ~QgramDiffer() = default;

  // amplist is the candidate list; difflist receives one distance per
  // candidate, so the two must carry the same length.
  //
  // The caller's underlying buffers are scratch sized to the whole pool,
  // and only their first few entries hold this seed's candidates. That is
  // why the length must come from the caller, who collected it, and not
  // from the buffer: reading a container's size() here was what commit
  // 799d763 did, and restoring the explicit length was half of the fix in
  // e517c04. A View carries the length the caller means -- the view *is*
  // the truncated list -- so the length stays explicit while travelling
  // with its data instead of alongside it.
  auto fast(uint64_t seed,
            View<uint64_t> amplist,
            Span<uint64_t> difflist) -> void;

  // The same, for candidates named by a run of amplicon pool entries rather
  // than by a list of ids the caller gathered first. The gather it replaces
  // was a pass over the whole pool per cluster, on the calling thread, to
  // copy four bytes out of each twenty-byte record -- 4 711 615 447
  // iterations and 1.5 s of a -d 2 run on 219k reads -- and the workers can
  // read those ids themselves while they are reading the q-gram vectors
  // anyway. The pool is read-only here: the distances still come back
  // through difflist, so nothing is written into the caller's records.
  auto fast_over_pool(uint64_t seed,
                      View<struct ampliconinfo_s> poollist,
                      Span<uint64_t> difflist) -> void;

private:
  // A thread's chunk names its candidates one of two ways: as a list of
  // amplicon ids (amplist, from fast()) or as a run of pool records to read
  // the ids out of (poollist, from fast_over_pool()). Exactly one is
  // non-empty, and worker() picks on that -- once per chunk, not once per
  // candidate.
  struct thread_info_s {
    uint64_t seed;
    View<uint64_t> amplist;                     // input, or empty
    View<struct ampliconinfo_s> poollist;       // input, or empty
    Span<uint64_t> difflist;  // output: the same extent, written by worker()
  };

  auto worker(uint64_t nth_thread) const noexcept -> void;

  // Hands each thread_info_s the [offset, chunk) slice of a list of
  // 'listlen' candidates that is its share, through 'assign', and runs them.
  // Shared by the two public entry points, which differ only in which member
  // of thread_info_s the slice lands in.
  template <typename Assign>
  auto distribute_and_run(uint64_t seed, uint64_t listlen, Assign assign) -> void;

  Qgram_store const          store_;        // owned (built in the ctor)
  Cpu_features const         cpu_features_;
  ThreadCount const          n_threads_;    // the configured -t value
  std::vector<thread_info_s> thread_info_v_;
  ThreadRunner               threads_;  // last: its lambda touches the members above
};

#endif  // SWARM_UTILS_QGRAM_H
