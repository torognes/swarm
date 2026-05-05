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

#include <cstdint>  // uint64_t
#include <vector>


// Blocked Bloom filter with precomputed bit patterns
// (Putze, Sanders, Singler 2009 -- see bloom.cc for the reference).
// Bit semantics are inverted from a textbook Bloom filter: a freshly
// zapped filter has all bits set to 1, set() clears the pattern's
// bits, and get() returns true (possibly-present) when all of the
// pattern's bits in the addressed word are zero.
class BloomFilter {
public:

  // bitmap_bytes is the requested bitmap size in bytes; it is rounded
  // up to at least one 64-bit word so bitmap_index() can compute a
  // valid position.
  BloomFilter(uint64_t bitmap_bytes,
              unsigned int shift,
              unsigned int n_hash_functions);

  // Reset all bitmap bits to 1 ("filter is empty"), preserving
  // allocated capacity and the precomputed patterns.
  auto zap() noexcept -> void;

  // Mark hash as a member of the set.
  auto set(uint64_t hash) noexcept -> void;

  // Test whether hash may be a member of the set. Returns true on
  // possible-membership, false on definite-non-membership.
  auto get(uint64_t hash) const noexcept -> bool;

private:

  auto bitmap_index(uint64_t hash) const noexcept -> uint64_t;
  auto bit_pattern(uint64_t hash) const noexcept -> uint64_t;
  auto generate_patterns() -> void;

  uint64_t size {0};            // bitmap length, in 64-bit words
  uint64_t pattern_shift {0};
  uint64_t pattern_count {0};
  uint64_t pattern_mask {0};
  uint64_t pattern_k {0};
  std::vector<uint64_t> bitmap;
  std::vector<uint64_t> patterns;
};
