/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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

#include "bloompat.h"
#include "utils/pseudo_rng.h"
#include <algorithm>  // std::max
#include <cstdint>  // uint64_t
#include <cstring>


void bloom_patterns_generate(struct bloom_s * bloom_filter)
{
  static constexpr unsigned int max_range {63};  // i & max_range = cap values to 63 max
  static constexpr unsigned int k {8};
  for(auto & pattern : bloom_filter->patterns)
    {
      pattern = 0;
      for(auto j = 0U; j < k; j++)
        {
          uint64_t onebit = 1ULL << (rand_64() & max_range);  // 0 <= shift <= 63
          while ((pattern & onebit) != 0) {
            onebit = 1ULL << (rand_64() & max_range);
          }
          pattern |= onebit;
        }
    }
}


void bloom_zap(struct bloom_s * bloom_filter)
{
  std::memset(bloom_filter->bitmap, UINT8_MAX, bloom_filter->size);
}


auto bloom_init(uint64_t size) -> struct bloom_s *
{
  // Size is in bytes for full bitmap, must be power of 2,
  // at least 8
  static constexpr uint64_t bytes_per_uint64 {8};
  size = std::max(size, bytes_per_uint64);

  auto * bloom_filter = new struct bloom_s;

  bloom_filter->size = size;

  bloom_filter->mask = (size >> 3U) - 1;

  bloom_filter->bitmap = new uint64_t[size];

  bloom_zap(bloom_filter);

  bloom_patterns_generate(bloom_filter);

  return bloom_filter;
}


void bloom_exit(struct bloom_s * bloom_filter)
{
  delete [] bloom_filter->bitmap;
  bloom_filter->bitmap = nullptr;
  delete bloom_filter;
  bloom_filter = nullptr;
}
