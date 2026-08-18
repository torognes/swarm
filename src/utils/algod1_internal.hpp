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
  extracted from algod1.cpp so that cohesive modules (network building,
  statistics, output) can live in their own translation units while
  still agreeing on the layout of the amplicon and swarm records.
*/

#include "../db.hpp"
#include "bloom.hpp"
#include "chain_range.hpp"  // Chain_range
#include "hashtable.hpp"
#include "span.hpp"  // Span, make_span
#include "view.hpp"  // View, make_view
#include <cassert>  // assert()
#include <cstddef>  // std::size_t
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
  // link_start is first so its 8-byte alignment does not pad the struct:
  // placed mid-struct it would grow ampinfo_s to 40 bytes instead of 32.
  uint64_t link_start {0U};            /* offset into network_v; total edges can exceed 2^32 */
  unsigned int swarmid {no_swarm};
  unsigned int parent {0U};
  unsigned int generation {0U};
  unsigned int next {no_swarm};        /* amp id of next amplicon in swarm */
  unsigned int graft_cand {no_swarm};  /* amp id of potential grafting parent (fastid.) */
  unsigned int link_count {0U};        /* per-amplicon, bounded by the variant count */
};


// The amplicons of one cluster are a chain through ampinfo_s::next, ending
// at no_swarm (see chain_range.hpp).
struct Next_cluster_member {
  static auto next(std::vector<struct ampinfo_s> const & ampinfo_v,
                   unsigned int const amp_id) -> unsigned int {
    return ampinfo_v[amp_id].next;
  }
};

using Cluster_chain = Chain_range<std::vector<struct ampinfo_s>, Next_cluster_member>;

// The cluster seeded by `seed`, seed included. Consumers that print a
// separator before every member except the first still compare against the
// seed themselves; the chain only supplies the order.
inline auto cluster_members(std::vector<struct ampinfo_s> const & ampinfo_v,
                            unsigned int const seed) -> Cluster_chain {
  return Cluster_chain{ampinfo_v, seed, no_swarm};
}

// The same cluster without its seed, for the consumers that report the seed
// separately and then walk the rest.
inline auto cluster_members_after_seed(std::vector<struct ampinfo_s> const & ampinfo_v,
                                       unsigned int const seed) -> Cluster_chain {
  return Cluster_chain{ampinfo_v, ampinfo_v[seed].next, no_swarm};
}


// One amplicon's slice of the flat network: its (link_start, link_count)
// pair resolved against network_v. Written once here rather than at each
// consumer, so that neither the offset nor the count can be applied to
// the wrong buffer, and so that the
// link_start + link_count <= network_v.size() bound is asserted (by
// subview/subspan) in every debug build instead of at one call site.
//
// Overloaded on the constness of network_v: the network writer sorts an
// amplicon's neighbours in place and needs a Span, while the cluster
// growth only reads them and takes a View.
inline auto neighbours_of(std::vector<unsigned int> const & network_v,
                          struct ampinfo_s const & amplicon) noexcept -> View<unsigned int> {
  return make_view(network_v).subview(static_cast<std::size_t>(amplicon.link_start),
                                     amplicon.link_count);
}

inline auto neighbours_of(std::vector<unsigned int> & network_v,
                          struct ampinfo_s const & amplicon) noexcept -> Span<unsigned int> {
  return make_span(network_v).subspan(static_cast<std::size_t>(amplicon.link_start),
                                     amplicon.link_count);
}

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
  uint64_t count {0};  /* running total of network edges; index into network_v */
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


inline auto hash_insert(Data const & data,
                        Hashtable & hash_table,
                        BloomFilter & bloom_a,
                        unsigned int const amp) -> void {
  /* find the first empty bucket */
  auto const hash = data.sequence_hash(amp);
  auto index = hash_table.getindex(hash);
  while (hash_table.is_occupied(index)) {
    index = hash_table.getnextindex(index);
  }

  hash_table.set_occupied(index);
  hash_table.set_value(index, hash);
  hash_table.set_data(index, amp);
  bloom_a.set(hash);
}

// max number of microvariants = 7 * len + 4
inline auto compute_microvariant_buffer_size(unsigned int const longest_sequence) noexcept -> unsigned int {
  static constexpr auto multiplier = 7U;
  static constexpr auto offset = 4U;
  // guard against unsigned int overflow of 7 * len + 4 + 1
  assert(longest_sequence <= (std::numeric_limits<unsigned int>::max() - offset - 1) / multiplier);
  return (multiplier * longest_sequence) + offset + 1;
}

#endif  // SWARM_UTILS_ALGOD1_INTERNAL_H
