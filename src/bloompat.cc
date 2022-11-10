/*
    SWARM

    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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

#include <cstdint>
#include <cstring>
#include "bloompat.h"
#include "pseudo_rng.h"


void bloom_patterns_generate(struct bloom_s * b)
{
  static constexpr unsigned int max_range {63};  // i & max_range = cap values to 63 max
  static constexpr unsigned int k {8};
  for(auto i = 0U; i < bloom_pattern_count; i++)  // range-loop '(unsigned long & i : b->patterns)' -> FAIL issue 123??
    {
      uint64_t pattern {0};
      for(auto j = 0U; j < k; j++)
        {
          uint64_t onebit {0};
          onebit = 1ULL << (rand_64() & max_range);  // 0 <= shift <= 63
          while ((pattern & onebit) != 0) {
            onebit = 1ULL << (rand_64() & max_range);
          }
          pattern |= onebit;
        }
      b->patterns[i] = pattern;
    }
}


void bloom_zap(struct bloom_s * b)
{
  memset(b->bitmap, UINT8_MAX, b->size);
}


auto bloom_init(uint64_t size) -> struct bloom_s *
{
  // Size is in bytes for full bitmap, must be power of 2,
  // at least 8
  static constexpr uint64_t bytes_per_uint64 {8};
  size = std::max(size, bytes_per_uint64);

  auto * b = new struct bloom_s;

  b->size = size;

  b->mask = (size >> 3) - 1;

  b->bitmap = new uint64_t[size];

  bloom_zap(b);

  bloom_patterns_generate(b);

  return b;
}


void bloom_exit(struct bloom_s * b)
{
  delete [] b->bitmap;
  b->bitmap = nullptr;
  delete b;
  b = nullptr;
}
