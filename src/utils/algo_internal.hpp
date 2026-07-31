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
  extracted from algo.cpp so that cohesive modules (currently the result
  writers in algo_output) can live in their own translation unit while
  still agreeing on the layout of the amplicon, swarm and workspace
  records.
*/

#include "span.hpp"  // Span, make_span
#include "view.hpp"  // View, make_view
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <vector>


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
  std::vector<uint64_t> alignlengths_v;
  std::vector<uint64_t> qgramamps_v;
  std::vector<uint64_t> qgramdiffs_v;
  std::vector<uint64_t> qgramindices_v;
  std::vector<uint64_t> hits;

  explicit Cluster_workspace(uint64_t const amplicons)
    : targetampliconids(amplicons),
      targetindices(amplicons),
      scores_v(amplicons),
      diffs_v(amplicons),
      alignlengths_v(amplicons),
      qgramamps_v(amplicons),
      qgramdiffs_v(amplicons),
      qgramindices_v(amplicons),
      hits(amplicons) {
  }

  // Every buffer above is allocated for the whole amplicon pool, while the
  // clustering code works on the first `count` entries of several of them
  // at once: the q-gram candidate list and its distances, then the target
  // list and its three result columns. The accessors below hand out those
  // windows, so a count is no longer paired with a buffer by hand at each
  // call site, and first() checks the window against the allocation.
  //
  // The read-only windows are const members; the ones the searches write
  // through are not, because a Span cannot come from a const container.
  auto qgram_candidates(std::size_t const count) const noexcept -> View<uint64_t> {
    return make_view(qgramamps_v).first(count);
  }
  auto qgram_diffs(std::size_t const count) noexcept -> Span<uint64_t> {
    return make_span(qgramdiffs_v).first(count);
  }
  auto targets(std::size_t const count) const noexcept -> View<uint64_t> {
    return make_view(targetampliconids).first(count);
  }
  auto scores(std::size_t const count) noexcept -> Span<uint64_t> {
    return make_span(scores_v).first(count);
  }
  auto diffs(std::size_t const count) noexcept -> Span<uint64_t> {
    return make_span(diffs_v).first(count);
  }
  auto alignlengths(std::size_t const count) noexcept -> Span<uint64_t> {
    return make_span(alignlengths_v).first(count);
  }
};

#endif  // SWARM_UTILS_ALGO_INTERNAL_H
