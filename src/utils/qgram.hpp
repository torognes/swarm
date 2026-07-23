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

#include "cpu_features.hpp"  // Cpu_features
#include "qgram_array.hpp"  // Qgram_store
#include "threads.hpp"  // ThreadRunner
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

  // listlen entries are read from the front of amplist; the same number
  // of entries are written to the front of difflist. amplist and difflist
  // are scratch buffers sized to the whole pool, so listlen (the number of
  // candidates actually collected) is passed explicitly rather than read
  // from amplist.size().
  auto fast(uint64_t seed,
            uint64_t listlen,
            std::vector<uint64_t> const & amplist,
            std::vector<uint64_t> & difflist) -> void;

private:
  struct thread_info_s {
    uint64_t seed;
    uint64_t listlen;
    uint64_t const * amplist;   // input: read-only inside worker()
    uint64_t * difflist;        // output: written by worker()
  };

  auto worker(uint64_t nth_thread) const noexcept -> void;

  Qgram_store const          store_;        // owned (built in the ctor)
  Cpu_features const         cpu_features_;
  std::vector<thread_info_s> thread_info_v_;
  ThreadRunner               threads_;  // last: its lambda touches the members above
};
