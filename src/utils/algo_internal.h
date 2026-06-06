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

#ifndef SWARM_UTILS_ALGO_INTERNAL_H
#define SWARM_UTILS_ALGO_INTERNAL_H


/*
  Shared records for the d > 1 clustering algorithm. These were
  extracted from algo.cc so that cohesive modules (currently the result
  writers in algo_output) can live in their own translation unit while
  still agreeing on the layout of the amplicon, swarm and workspace
  records.
*/

#include <cstdint>  // int64_t, uint64_t
#include <vector>

struct Parameters;    // defined in swarm.h
class Data;           // defined in db.h
class QgramDiffer;    // defined in utils/qgram.h
struct Search_state;  // defined in scan.h
class ThreadRunner;   // defined in utils/threads.h


struct ampliconinfo_s {
  unsigned int ampliconid {0};
  unsigned int diffestimate {0}; /* lower bound estimate of dist from initial seed */
  unsigned int swarmid {0};
  unsigned int generation {0};
  unsigned int radius {0}; /* actual diff from initial seed */
};

struct swarminfo_t {
  uint64_t mass {0};
  unsigned int seed {0};
  int dummy {0}; /* alignment padding only */
};


struct Cluster_state {
  uint64_t swarmsize {1};         // a cluster cannot be empty
  uint64_t amplicons_copies {0};  // total abundance of the cluster
  uint64_t singletons {0};
  uint64_t hitcount {0};
  uint64_t maxradius {0};
  uint64_t maxgen {1};            // a cluster can't contain less than 1 generation
};


// Cursors into the amplicon pool, shared across all clusters. Bundled
// into a struct so they are never passed as two swappable uint64_t
// arguments to the clustering helpers.
struct Pool_cursor {
  uint64_t seeded {0};   // amplicons whose subseed expansion is complete
  uint64_t swarmed {0};  // amplicons already assigned to a swarm
};


struct Cluster_workspace {
  std::vector<uint64_t> targetampliconids;
  std::vector<uint64_t> targetindices;
  std::vector<uint64_t> scores_v;
  std::vector<uint64_t> diffs_v;
  std::vector<uint64_t> alignlengths;
  std::vector<uint64_t> qgramamps_v;
  std::vector<uint64_t> qgramdiffs_v;
  std::vector<uint64_t> qgramindices_v;
  std::vector<uint64_t> hits;

  explicit Cluster_workspace(uint64_t const amplicons)
    : targetampliconids(amplicons),
      targetindices(amplicons),
      scores_v(amplicons),
      diffs_v(amplicons),
      alignlengths(amplicons),
      qgramdiffs_v(amplicons),
      hits(amplicons) {
    qgramamps_v.reserve(amplicons);
    qgramindices_v.reserve(amplicons);
  }
};


// Read-only bundle of the search machinery shared by every cluster:
// the qgram pre-filter, the alignment search engine and its worker
// threads. Holds references/pointers only, so it is cheap to pass and
// owns nothing. (search_state and qgram_differ are mutated through
// their references; the surrounding const only fixes the bindings.)
struct Search_context {
  struct Parameters const & parameters;
  Data const & data;
  QgramDiffer & qgram_differ;
  struct Search_state & search_state;
  ThreadRunner * search_threads;
  int bits;
};

#endif  // SWARM_UTILS_ALGO_INTERNAL_H
