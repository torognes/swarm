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

#ifndef SWARM_UTILS_ALGOD1_INTERNAL_H
#define SWARM_UTILS_ALGOD1_INTERNAL_H


/*
  Shared types, constants and helpers for the d=1 algorithm. These were
  extracted from algod1.cc so that cohesive modules (network building,
  statistics, output) can live in their own translation units while
  still agreeing on the layout of the amplicon and swarm records.
*/

#include <cstdint>  // int64_t, uint64_t
#include <limits>  // unsigned int max
#include <mutex>  // std::mutex
#include <vector>


constexpr unsigned int one_kilobyte {1U << 10U};  // 1,024 bytes
constexpr unsigned int one_megabyte {one_kilobyte * one_kilobyte};
constexpr unsigned int no_swarm {std::numeric_limits<unsigned int>::max()};


/* Information about each amplicon */

struct ampinfo_s
{
  unsigned int swarmid {no_swarm};
  unsigned int parent {0U};
  unsigned int generation {0U};
  unsigned int next {no_swarm};        /* amp id of next amplicon in swarm */
  unsigned int graft_cand {no_swarm};  /* amp id of potential grafting parent (fastid.) */
  unsigned int link_start {0U};
  unsigned int link_count {0U};
};

/* Information about each swarm (cluster) */

struct swarminfo_s
{
  uint64_t mass {0}; /* the sum of abundances of amplicons in this swarm */
  uint64_t sumlen {0}; /* sum of length of amplicons in swarm */
  unsigned int seed {0}; /* amplicon id of the initial seed of this swarm */
  unsigned int last {0}; /* amplicon id of the last seed in this swarm */
  unsigned int size {0}; /* total number of amplicons in this swarm */
  unsigned int singletons {0}; /* number of amplicons with abundance 1 */
  unsigned int maxgen {0}; /* the generation of the amplicon farthest from seed */
  bool attached {false}; /* this is a small swarm attached to a large (fastidious) */
  char dummy_1 = '\0'; /* alignment padding only */
  char dummy_2 = '\0'; /* alignment padding only */
  char dummy_3 = '\0'; /* alignment padding only */
};  // total of 40 bytes (five 64-bit machine words)

/* overall statistics, accumulated across all swarms */
struct Overall_stats
{
  uint64_t swarmcount_adjusted {0};
  unsigned int maxgen {0};
  unsigned int largest {0};
};

struct Network_state
{
  std::mutex mutex;
  unsigned int amp {0};
  unsigned int count {0};
  std::vector<unsigned int> network_v;
};

struct Cluster_stats
{
  uint64_t small_clusters {0};
  uint64_t large_clusters {0};
  uint64_t amplicons_in_small_clusters {0};
  uint64_t amplicons_in_large_clusters {0};
  uint64_t nucleotides_in_small_clusters {0};
};

struct Bloom_geometry
{
  uint64_t n_bytes {0};
  unsigned int n_hash_functions {0};
};

/* Bloom filter shape used for the per-amplicon hashtable + bloom_a
   in both the d=1 phase and the fastidious phase. */
constexpr unsigned int amplicon_pattern_shift {10};
constexpr unsigned int amplicon_n_hash_functions {8};

#endif  // SWARM_UTILS_ALGOD1_INTERNAL_H
