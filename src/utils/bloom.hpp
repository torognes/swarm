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


#ifndef SWARM_UTILS_BLOOM_H
#define SWARM_UTILS_BLOOM_H

#include <algorithm>  // std::max
#include <cassert>
#include <cstdint>  // uint64_t
#include <limits>  // std::numeric_limits
#include <vector>


namespace bloom_detail {

  // bitmap is stored as an array of 64-bit words; this is the size of
  // one such word (in bytes), used both to lower-bound the requested
  // size and to convert bytes -> words via a right shift.
  constexpr uint64_t bytes_per_word {8};
  static_assert(bytes_per_word == sizeof(uint64_t),
                "bytes_per_word must match sizeof(uint64_t)");

  // Give each pattern n_hash_functions distinct bits, drawn at random.
  //
  // Not a member and not a template: the pattern *shift* decides how many
  // patterns there are, but nothing here depends on it, so this loop is
  // compiled once for the whole program instead of once per instantiation
  // -- and <random> stays out of a header. Defined in bloom.cpp.
  //
  // Not marked noexcept: rand_64.operator() (std::mt19937_64) is not
  // formally noexcept in the standard, even if it does not throw in
  // practice. Called only from the constructor, which is itself non-
  // noexcept (vector resize/construction can throw bad_alloc), so the
  // distinction is academic.
  auto generate_patterns(std::vector<uint64_t> & patterns,
                         uint64_t pattern_k) -> void;

}  // namespace bloom_detail


// Blocked Bloom filter with precomputed bit patterns
// (Putze, Sanders, Singler 2009 -- see bloom.cpp for the reference).
// Bit semantics are inverted from a textbook Bloom filter: a freshly
// constructed filter has all bits set to 1, set() clears the pattern's
// bits, and get() returns true (possibly-present) when all of the
// pattern's bits in the addressed word are zero.
//
// pattern_shift is a template parameter rather than a constructor one:
// every filter swarm builds names it with a constant (amplicon_pattern_shift
// and fastidious_pattern_shift, both in algod1_internal.hpp), so the
// runtime parameter only ever carried a compile-time fact. Lifting it into
// the type turns the shift in bitmap_index() and the mask in bit_pattern()
// from loaded members into immediates, and removes four of the five scalar
// members: the shift itself, the count and mask derived from it, and
// pattern_k, which only the pattern generation above ever read. Both are on
// the probe path, which runs twice per Bloom lookup in the fastidious
// passes.
//
// It also gives the two filters distinct types. Six functions in
// algod1_fastidious.cpp take bloom_a and bloom_f as adjacent parameters,
// which were the same type and so could be passed in either order; now
// they cannot.
//
// n_hash_functions stays a constructor parameter: the fastidious filter
// derives its k from --bloom-bits (see compute_bloom_geometry), and k is
// read only by the cold pattern generation, so there would be nothing to
// fold.
template <unsigned int pattern_shift>
class BloomFilter {
public:

  // bitmap_bytes is the requested bitmap size in bytes; it is rounded
  // up to at least one 64-bit word so bitmap_index() can compute a
  // valid position.
  //
  // Non-noexcept: the two vector constructions from (count, value) can
  // throw std::bad_alloc.
  BloomFilter(uint64_t const bitmap_bytes,
              unsigned int const n_hash_functions)
    : size{std::max(bitmap_bytes, bloom_detail::bytes_per_word) >> 3U}
    , bitmap(size, std::numeric_limits<uint64_t>::max())
    , patterns(pattern_count) {
    bloom_detail::generate_patterns(patterns, n_hash_functions);
  }

  // Mark hash as a member of the set.
  auto set(uint64_t const hash) noexcept -> void {
    bitmap[bitmap_index(hash)] &= compl bit_pattern(hash);
  }

  // Test whether hash may be a member of the set. Returns true on
  // possible-membership, false on definite-non-membership.
  auto get(uint64_t const hash) const noexcept -> bool {
    return (bitmap[bitmap_index(hash)] & bit_pattern(hash)) == 0U;
  }

private:

  // Derived from the shift, so constants rather than members.
  static constexpr uint64_t pattern_count {uint64_t{1} << pattern_shift};
  static constexpr uint64_t pattern_mask {pattern_count - 1};

  // Refactoring: the modulo below is on the hot path (called twice per
  // Bloom filter probe in algod1.cpp) and is markedly slower than a
  // bitwise AND. The previous bloompat code used `& mask` because it
  // required `size` to be a power of 2; BloomFilter accepts arbitrary
  // sizes, so it must use `%`. To restore the fast path, constrain
  // `size` to be a power of 2 (round up or down in the constructor or
  // in the caller), store `size - 1` as a mask, and replace `% size`
  // with `& mask`. The amplicon filter already receives a power-of-2
  // size from compute_hashtable_size(); the fastidious filter does
  // not, and would need its caller in algod1.cpp to choose a rounding
  // policy compatible with the --ceiling / --bloom-bits memory budget.
  //
  auto bitmap_index(uint64_t const hash) const noexcept -> uint64_t {
    auto const position = (hash >> pattern_shift) % size;
    assert(position < bitmap.size());
    return position;
  }

  auto bit_pattern(uint64_t const hash) const noexcept -> uint64_t {
    auto const position = hash & pattern_mask;
    assert(position < patterns.size());
    return patterns[position];
  }

  uint64_t size {0};            // bitmap length, in 64-bit words
  std::vector<uint64_t> bitmap;
  std::vector<uint64_t> patterns;
};

#endif  // SWARM_UTILS_BLOOM_H
