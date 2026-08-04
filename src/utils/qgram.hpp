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
#include "qgram_array.hpp"  // Qgram_store
#include "span.hpp"  // Span<uint64_t>
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

private:
  struct thread_info_s {
    uint64_t seed;
    View<uint64_t> amplist;   // input: this thread's chunk, read-only in worker()
    Span<uint64_t> difflist;  // output: the same extent, written by worker()
  };

  auto worker(uint64_t nth_thread) const noexcept -> void;

  Qgram_store const          store_;        // owned (built in the ctor)
  Cpu_features const         cpu_features_;
  std::vector<thread_info_s> thread_info_v_;
  ThreadRunner               threads_;  // last: its lambda touches the members above
};

#endif  // SWARM_UTILS_QGRAM_H
