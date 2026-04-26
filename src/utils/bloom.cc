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

/*
  Blocked bloom filter with precomputed bit patterns
  as described in

  Putze F, Sanders P, Singler J (2009)
  Cache-, Hash- and Space-Efficient Bloom Filters
  Journal of Experimental Algorithmics, 14, 4
  https://doi.org/10.1145/1498698.1594230
*/

#include "bloom.h"
#include "pseudo_rng.h"
#include <algorithm>  // std::fill, std::max
#include <cassert>
#include <cstdint>  // uint64_t
#include <limits>


// Constructor is non-noexcept: the two vector resizes / construction
// from (count, value) can throw std::bad_alloc.
BloomFilter::BloomFilter(const uint64_t bitmap_bytes,
                         const unsigned int shift,
                         const unsigned int n_hash_functions)
  : size{std::max(bitmap_bytes, uint64_t{8}) >> 3U}
  , pattern_shift{shift}
  , pattern_count{uint64_t{1} << shift}
  , pattern_mask{pattern_count - 1}
  , pattern_k{n_hash_functions}
  , bitmap(size, std::numeric_limits<uint64_t>::max())
  , patterns(pattern_count)
{
  generate_patterns();
}


// Refactoring: the modulo below is on the hot path (called twice per
// Bloom filter probe in algod1.cc) and is markedly slower than a
// bitwise AND. The previous bloompat code used `& mask` because it
// required `size` to be a power of 2; BloomFilter accepts arbitrary
// sizes, so it must use `%`. To restore the fast path, constrain
// `size` to be a power of 2 (round up or down in the constructor or
// in the caller), store `size - 1` as a mask, and replace `% size`
// with `& mask`. The amplicon filter already receives a power-of-2
// size from compute_hashtable_size(); the fastidious filter does
// not, and would need its caller in algod1.cc to choose a rounding
// policy compatible with the --ceiling / --bloom-bits memory budget.
//
// noexcept: arithmetic plus an unchecked vector subscript is implicit
// here only via the assert; the body never allocates or throws.
auto BloomFilter::bitmap_index(const uint64_t hash) const noexcept -> uint64_t
{
  auto const position = (hash >> pattern_shift) % size;
  assert(position < bitmap.size());
  return position;
}


// noexcept: arithmetic plus a bounds-checked vector subscript;
// operator[] does not throw and the assert guards out-of-range access
// in debug builds.
auto BloomFilter::pat(const uint64_t hash) const noexcept -> uint64_t
{
  auto const position = hash & pattern_mask;
  assert(position < patterns.size());
  return patterns[position];
}


// noexcept: only calls noexcept helpers and performs vector subscript
// (does not throw) plus bitwise arithmetic on built-ins.
auto BloomFilter::set(const uint64_t hash) noexcept -> void
{
  bitmap[bitmap_index(hash)] &= compl pat(hash);
}


// noexcept: same reasoning as set().
auto BloomFilter::get(const uint64_t hash) const noexcept -> bool
{
  return (bitmap[bitmap_index(hash)] & pat(hash)) == 0U;
}


// noexcept: std::fill on a std::vector with already-allocated storage
// performs no allocation and only invokes uint64_t assignment, which
// cannot throw.
auto BloomFilter::zap() noexcept -> void
{
  std::fill(bitmap.begin(), bitmap.end(), std::numeric_limits<uint64_t>::max());
}


// Not marked noexcept: rand_64.operator() (std::mt19937_64) is not
// formally noexcept in the standard, even if it does not throw in
// practice. Called only from the constructor, which is itself non-
// noexcept (vector resize/construction can throw bad_alloc), so the
// distinction is academic.
auto BloomFilter::generate_patterns() -> void
{
  static constexpr auto max_range = 63U;  // i & max_range = cap values to 63 max
  for(auto & pattern : patterns)
    {
      pattern = 0;
      for(auto j = 0U; j < pattern_k; ++j)
        {
          uint64_t onebit = 1ULL << (rand_64() & max_range);  // 0 <= shift <= 63
          while ((pattern & onebit) != 0U) {
            onebit = 1ULL << (rand_64() & max_range);
          }
          pattern |= onebit;
        }
    }
}
