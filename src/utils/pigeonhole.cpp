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

#include "pigeonhole.hpp"
#include "../db.hpp"  // Data, Sequence, nucleotide_at
#include "../swarm.hpp"  // struct Parameters
#include "memory_budget.hpp"  // require_ram
#include "view.hpp"  // make_view
#include <algorithm>  // std::min, std::sort, std::fill, std::is_sorted
#include <cassert>
#include <cstdint>  // uint64_t, int64_t
#include <iterator>  // std::next
#include <utility>  // std::pair
#include <vector>


namespace {

  struct Segment_bounds {
    uint64_t start;
    uint64_t length;
  };

  // Segment number 'segment' (0-based) of a sequence of 'length'
  // nucleotides cut into 'n_segments' near-equal parts: the first
  // (length % n_segments) segments are one nucleotide longer than the
  // rest. Build and query must agree on this rule, which is why it lives
  // in one place.
  auto segment_bounds(uint64_t const length,
                      uint64_t const segment,
                      uint64_t const n_segments) -> Segment_bounds {
    auto const base = length / n_segments;
    auto const n_longer = length % n_segments;
    auto const start = (segment * base) + std::min(segment, n_longer);
    auto const seg_length = base + (segment < n_longer ? 1 : 0);
    return {start, seg_length};
  }


  // The segment content is hashed with a polynomial in this multiplier
  // over the 2-bit nucleotide codes, mod 2^64: unlike the FNV-1a it
  // replaced, a polynomial rolls -- shifting the window one position
  // costs one subtract-multiply-add instead of a pass over the segment,
  // which is what search()'s shift loop needs. The multiplier is
  // Knuth's 64-bit LCG constant; anything odd and well mixed does, since
  // a splitmix64 finalizer spreads the result afterwards and collisions
  // are safe by design (see the class comment).
  constexpr uint64_t poly_multiplier {6364136223846793005ULL};

  // The content polynomial of sequence[start, start + seg_length):
  //   code[start] * B^(seg_length - 1) + ... + code[start + seg_length - 1]
  auto segment_poly(Sequence const & sequence,
                    uint64_t const start,
                    uint64_t const seg_length) -> uint64_t {
    uint64_t poly {0};
    auto const past_last = start + seg_length;
    for (auto position = start; position < past_last; ++position) {
      poly = (poly * poly_multiplier) + nucleotide_at(sequence, position);
    }
    return poly;
  }

  // The bucket key of one (target length, segment number, content)
  // triple: the coordinates folded into the content polynomial, then a
  // splitmix64 finalizer so that every input bit reaches the low bits
  // the bucket mask keeps.
  auto finalize_key(uint64_t const poly,
                    uint64_t const target_length,
                    uint64_t const segment) -> uint64_t {
    // splitmix64 finalizer constants (Steele, Lea & Flood 2014)
    static constexpr uint64_t mix_multiplier_1 {0xbf58476d1ce4e5b9ULL};
    static constexpr uint64_t mix_multiplier_2 {0x94d049bb133111ebULL};
    static constexpr unsigned int mix_shift_1 {30};
    static constexpr unsigned int mix_shift_2 {27};
    static constexpr unsigned int mix_shift_3 {31};
    // the segment number occupies the low bits, the length the rest;
    // segments are numbered 0..d with d < 256
    static constexpr unsigned int segment_number_bits {8};

    auto key = poly ^ (target_length << segment_number_bits) ^ segment;
    key ^= key >> mix_shift_1;
    key *= mix_multiplier_1;
    key ^= key >> mix_shift_2;
    key *= mix_multiplier_2;
    key ^= key >> mix_shift_3;
    return key;
  }

}  // namespace


PigeonholeIndex::PigeonholeIndex(struct Parameters const & parameters,
                                 Data const & data)
  : data_(data),
    n_differences_(parameters.opt_differences)
{
  auto const n_sequences = uint64_t{data.sequence_count()};
  auto const n_segments = n_differences_ + 1;
  auto const n_entries = n_sequences * n_segments;

  // the smallest power of two giving at most one entry per bucket on
  // average; unrelated triples sharing a bucket only add candidates
  uint64_t n_buckets {1};
  while (n_buckets < n_entries) {
    n_buckets += n_buckets;
  }
  bucket_mask_ = n_buckets - 1;

  // entries (4 bytes each), offsets and the build-time cursors (8 bytes
  // per bucket each, with fewer than two buckets per entry), and the
  // per-sequence query stamps (4 bytes, at most one sequence per entry):
  // under 40 bytes per entry all told. The product cannot overflow:
  // n_entries is at most 2^32 sequences times 256 segments, well under
  // 2^64 / 40.
  static constexpr uint64_t bytes_per_entry {40};
  require_ram(bytes_per_entry * n_entries, 1, "the pigeonhole index");

  offsets_.assign(n_buckets + 1, 0);
  entries_.resize(n_entries);
  length_present_.assign(data.longest_sequence() + n_differences_ + 1, false);
  seen_stamp_.assign(n_sequences, 0);

  // B^i for i up to the longest segment length, so search() can roll its
  // hash window: dropping the leading code costs one multiply by B^(len-1)
  multiplier_powers_.resize(uint64_t{data.longest_sequence()} + 1);
  multiplier_powers_[0] = 1;
  for (auto power = 1ULL; power < multiplier_powers_.size(); ++power) {
    multiplier_powers_[power] = multiplier_powers_[power - 1] * poly_multiplier;
  }

  // count the population of each bucket...
  for (auto seqno = 0ULL; seqno < n_sequences; ++seqno) {
    auto const sequence = data.sequence_view(seqno);
    length_present_[sequence.length] = true;
    for (auto segment = 0ULL; segment < n_segments; ++segment) {
      auto const bounds = segment_bounds(sequence.length, segment, n_segments);
      auto const key = finalize_key(segment_poly(sequence, bounds.start,
                                                 bounds.length),
                                    sequence.length, segment);
      ++offsets_[(key & bucket_mask_) + 1];
    }
  }

  // ...turn the counts into bucket boundaries...
  for (auto bucket = 1ULL; bucket <= n_buckets; ++bucket) {
    offsets_[bucket] += offsets_[bucket - 1];
  }
  assert(offsets_[n_buckets] == n_entries);

  // ...and place the amplicon ids. Ids are visited in ascending order, so
  // each bucket's ids come out ascending: search() counts on that to hand
  // back candidates in pool order after one merge-free sort.
  std::vector<uint64_t> cursors(offsets_.cbegin(),
                                std::next(offsets_.cbegin(),
                                          static_cast<std::ptrdiff_t>(n_buckets)));
  for (auto seqno = 0ULL; seqno < n_sequences; ++seqno) {
    auto const sequence = data.sequence_view(seqno);
    for (auto segment = 0ULL; segment < n_segments; ++segment) {
      auto const bounds = segment_bounds(sequence.length, segment, n_segments);
      auto const key = finalize_key(segment_poly(sequence, bounds.start,
                                                 bounds.length),
                                    sequence.length, segment);
      entries_[cursors[key & bucket_mask_]] = static_cast<unsigned int>(seqno);
      ++cursors[key & bucket_mask_];
    }
  }
}


auto PigeonholeIndex::search(uint64_t const query) -> Search_result {
  // Refuse to gather more candidates than this, and report the query
  // 'heavy' instead: the caller's fallback is the full pool scan, whose
  // per-candidate cost (a 128-byte q-gram comparison, parallelized) is
  // several times below this path's gather-sort-deduplicate cost per
  // candidate, so an index win requires far fewer candidates than the
  // pool holds. The value is measured on the d = 2 / d = 3 profiling
  // dataset (see the cap sweep in the d > 1 notes): it keeps 72-87 % of
  // d = 2 queries and 40-49 % of d = 3 queries on the index at a few
  // hundred candidates each, while the skewed rest -- queries touching a
  // conserved-region bucket -- pay one probe pass (~30-80 constant-time
  // size lookups) and fall back.
  static constexpr uint64_t max_gather_volume {4096};

  auto const & data = data_.get();
  auto const sequence = data.sequence_view(query);
  auto const query_length = static_cast<int64_t>(sequence.length);
  auto const n_segments = n_differences_ + 1;
  auto const differences = static_cast<int64_t>(n_differences_);

  probes_v_.clear();
  uint64_t volume {0};

  for (auto target_length = query_length - differences;
       target_length <= query_length + differences;
       ++target_length) {
    if (target_length < 0) { continue; }
    assert(static_cast<uint64_t>(target_length) < length_present_.size());
    if (not length_present_[static_cast<uint64_t>(target_length)]) { continue; }
    auto const delta = query_length - target_length;

    for (auto segment = 0ULL; segment < n_segments; ++segment) {
      auto const bounds =
        segment_bounds(static_cast<uint64_t>(target_length), segment, n_segments);
      auto const start = static_cast<int64_t>(bounds.start);
      auto const seg_length = static_cast<int64_t>(bounds.length);

      // A target within distance d has at least one segment untouched by
      // any edit, so that segment occurs verbatim in the query, displaced
      // by the net indel count of the edits before it. Those edits split
      // around the segment: e_pre before it bounds the displacement,
      // |shift| <= e_pre; e_post after it must absorb the rest of the
      // length difference, |delta - shift| <= e_post; and e_pre + e_post
      // <= d. So |shift| + |delta - shift| <= d, which is the contiguous
      // window (delta - d) / 2 <= shift <= (delta + d) / 2 -- with
      // |delta| <= d the left bound's numerator is never positive and the
      // right one's never negative, so truncating division lands on the
      // ceiling and the floor. (This is the length-aware window only --
      // the tighter multi-match-aware window of the PassJoin paper trades
      // a subtle proof for a constant factor, and a mistake there would
      // lose true neighbours; this one cannot.)
      auto const shift_low = std::max((delta - differences) / 2, -start);
      auto const shift_high = std::min((delta + differences) / 2,
                                       query_length - seg_length - start);
      if (shift_low > shift_high) { continue; }

      // The window is contiguous, so the content hash rolls: one full
      // pass for the first position, then one subtract-multiply-add per
      // shift. An empty segment hashes the same everywhere and is probed
      // once.
      auto position = static_cast<uint64_t>(start + shift_low);
      auto const last_position = static_cast<uint64_t>(start + shift_high);
      auto poly = segment_poly(sequence, position, bounds.length);
      while (true) {
        auto const key = finalize_key(poly,
                                      static_cast<uint64_t>(target_length),
                                      segment);
        auto const bucket = key & bucket_mask_;
        auto const begin = offsets_[bucket];
        auto const end = offsets_[bucket + 1];
        if (begin != end) {
          probes_v_.emplace_back(begin, end);
          volume += end - begin;
        }
        if ((position >= last_position) or (bounds.length == 0)) { break; }
        poly -= static_cast<uint64_t>(nucleotide_at(sequence, position))
              * multiplier_powers_[bounds.length - 1];
        poly = (poly * poly_multiplier)
             + nucleotide_at(sequence, position + bounds.length);
        ++position;
      }
    }
  }

  if (volume > max_gather_volume) {
    return {true, View<unsigned int>{}};
  }

  // a fresh stamp marks this query's round; when the counter wraps, the
  // whole array is cleared so no stamp left by an old round can alias it
  ++query_stamp_;
  if (query_stamp_ == 0) {
    std::fill(seen_stamp_.begin(), seen_stamp_.end(), 0U);
    ++query_stamp_;
  }

  candidates_v_.clear();
  candidates_v_.reserve(volume);
  for (auto const & probe : probes_v_) {
    for (auto index = probe.first; index < probe.second; ++index) {
      auto const candidate = entries_[index];
      if (seen_stamp_[candidate] == query_stamp_) { continue; }
      seen_stamp_[candidate] = query_stamp_;
      candidates_v_.push_back(candidate);
    }
  }
  // the stamps already deduplicated, so the sort -- pool order is part of
  // the contract -- runs on the unique candidates only
  std::sort(candidates_v_.begin(), candidates_v_.end());
  return {false, make_view(candidates_v_)};
}
