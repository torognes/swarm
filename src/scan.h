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

#include "utils/queryinfo.h"
#include "utils/threads.h"
#include <cstdint>  // int64_t, uint64_t
#include <mutex>
#include <vector>


struct Search_data;  // defined in utils/search_data.h
struct Parameters;  // defined in swarm.h
class  Data;        // defined in db.h

struct Search_state
{
  std::mutex scan_mutex;
  struct Search_data * search_data {nullptr};
  struct queryinfo query {0, 0, nullptr};
  uint64_t master_next {0};
  uint64_t master_length {0};
  uint64_t remainingchunks {0};
  uint64_t * master_targets {nullptr};
  uint64_t * master_scores {nullptr};
  uint64_t * master_diffs {nullptr};
  uint64_t * master_alignlengths {nullptr};
  int master_bits {0};
};


auto search_do(struct Parameters const & parameters,
               Data const & data,
               struct Search_state & state,
               uint64_t query_no,
               uint64_t listlength,
               uint64_t * targets,
               uint64_t * scores,
               uint64_t * diffs,
               uint64_t * alignlengths,
               int bits,
               ThreadRunner * search_threads) -> void;
auto search_begin(struct Parameters const & parameters,
                  Data const & data,
                  struct Search_state & state,
                  std::vector<struct Search_data> & search_data_v) -> void;
auto search_end(struct Search_state & state) -> void;
auto search_worker_core(struct Parameters const & parameters,
                        Data const & data,
                        uint64_t thread_id, struct Search_state & state) -> void;
