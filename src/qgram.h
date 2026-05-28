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

#include "utils/cpu_features.h"  // Cpu_features
#include "utils/qgram_array.h"  // Qgram_store
#include "utils/qgram_threadinfo.h"  // thread_info_s
#include "utils/threads.h"  // ThreadRunner
#include <cstdint>  // uint64_t
#include <vector>


struct Parameters;  // defined in swarm.h
class Data;          // defined in db.h

auto findqgrams(char const * seq, uint64_t seqlen,
                unsigned char * qgramvector) -> void;
auto build_qgram_store(struct Parameters const & parameters,
                       Data const & data) -> Qgram_store;
auto qgram_diff_fast(struct Parameters const & parameters,
                     Qgram_store const & store,
                     uint64_t seed,
                     uint64_t listlen,
                     uint64_t * amplist,
                     uint64_t * difflist,
                     std::vector<struct thread_info_s>& thread_info_v) -> void;
auto qgram_diff_init(struct Parameters const & parameters,
                     Qgram_store const & store,
                     std::vector<struct thread_info_s>& thread_info_v) -> void;
auto qgram_diff_done() -> void;


// RAII wrapper around the per-seed qgram-distance dispatch: owns the
// worker ThreadRunner, the per-thread thread_info_s scratch vector,
// and the cached Cpu_features. The destructor joins the worker
// threads, replacing the qgram_diff_init / qgram_diff_done pair.
class QgramDiffer {
public:
  QgramDiffer(struct Parameters const & parameters,
              Qgram_store const & store);

  // Non-copyable, non-movable: the ThreadRunner's lambda captures
  // `this`, so the object must keep a stable address.
  QgramDiffer(QgramDiffer const &) = delete;
  QgramDiffer(QgramDiffer &&) = delete;
  auto operator=(QgramDiffer const &) -> QgramDiffer & = delete;
  auto operator=(QgramDiffer &&) -> QgramDiffer & = delete;

  auto fast(uint64_t seed,
            uint64_t listlen,
            uint64_t * amplist,
            uint64_t * difflist) -> void;

private:
  Qgram_store const &              store_;
  Cpu_features const               cpu_features_;
  std::vector<struct thread_info_s> thread_info_v_;
  ThreadRunner                     threads_;  // last: its lambda touches the members above
};
