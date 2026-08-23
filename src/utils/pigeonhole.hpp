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

#ifndef SWARM_UTILS_PIGEONHOLE_H
#define SWARM_UTILS_PIGEONHOLE_H


#include "view.hpp"  // View<unsigned int>
#include <cstdint>  // uint64_t, std::uint32_t
#include <functional>  // std::reference_wrapper
#include <utility>  // std::pair
#include <vector>


struct Parameters;  // defined in swarm.hpp
class Data;         // defined in db.hpp


// An exact candidate index for the d > 1 clustering path, by pigeonhole /
// segment filtering (the PassJoin family: Li, Deng & Feng, ICDE 2011).
//
// Build: every sequence is cut into d + 1 segments of near-equal length,
// and each (sequence length, segment number, segment content) triple is
// indexed. Two sequences within edit distance d must agree exactly on at
// least one segment -- d edits cannot touch d + 1 segments -- so a query
// that probes the right triples cannot miss a true neighbour. Segment
// content enters the index as a 64-bit hash: a collision merges buckets,
// which can only add candidates, never hide one, and every candidate is
// verified downstream (q-gram bound, then alignment) exactly as
// scan-collected candidates are. The index is therefore exact in the only
// direction that matters -- no false negatives -- while staying free to
// return false positives.
//
// Query: search() enumerates, for each target length within d of the
// query's and each segment of that length's partition, the query
// substrings that segment could match, and gathers the probed buckets.
// The candidate volume is known from the bucket sizes before anything is
// gathered; a query whose buckets hold more than a full pool scan is
// worth is refused ('heavy') so the caller can run the scan instead.
// That bounds the damage of skew -- with segments this short, one
// conserved region can put a large fraction of the dataset in one bucket
// (measured: 6.5 % of a 219k-read V9 dataset behind a single d = 2 key).
class PigeonholeIndex {
public:
  // The largest d the index is built for. Query cost grows with d on
  // every axis -- (2d+1) target lengths x (d+1) segments x up to (d+1)
  // shifts, on segments that get shorter and buckets that get bigger --
  // and the measured gather volume goes from winning at d = 2 (0.7 G
  // candidates versus 6.2 G q-gram comparisons plus a 19.3 G-step serial
  // walk) to marginal at d = 3 (2.8 G) to hopeless at d = 10 (40 G, all
  // queries over any workable cap). d = 2 and d = 3 are also what the
  // downstream pipelines ask for; larger d stays on the scan path.
  static constexpr uint64_t max_indexed_differences {3};

  PigeonholeIndex(struct Parameters const & parameters, Data const & data);

  struct Search_result {
    // gathering was refused: the probed buckets hold more candidates
    // than the full scan this query would replace -- run the scan
    bool heavy;
    View<unsigned int> candidates;  // ascending amplicon ids, deduplicated
  };

  // Every amplicon within edit distance opt_differences of amplicon
  // 'query' is in candidates (the query itself included -- it matches its
  // own segments), along with false positives for the verification stages
  // to reject. The view aliases internal scratch: it is valid until the
  // next search() call, and search() must not be called concurrently.
  auto search(uint64_t query) -> Search_result;

private:
  std::reference_wrapper<Data const> data_;
  uint64_t n_differences_;
  uint64_t bucket_mask_ {0};
  std::vector<uint64_t> offsets_;      // bucket b: entries_[offsets_[b], offsets_[b+1])
  std::vector<unsigned int> entries_;  // amplicon ids, ascending within a bucket
  std::vector<bool> length_present_;   // which sequence lengths occur
  std::vector<std::pair<uint64_t, uint64_t>> probes_v_;  // search() scratch
  std::vector<unsigned int> candidates_v_;               // search() scratch and result
  // per-amplicon stamp of the last query that gathered it, so duplicates
  // -- one candidate matching on several segments -- are dropped at
  // gather time for one array load each, instead of surviving into the
  // sort and being erased by std::unique afterwards
  std::vector<std::uint32_t> seen_stamp_;
  std::uint32_t query_stamp_ {0};
};

#endif  // SWARM_UTILS_PIGEONHOLE_H
