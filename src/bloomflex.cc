/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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

#include "bloomflex.h"
#include "util.h"
#include "utils/pseudo_rng.h"
#include <cstdint>  // uint64_t
#include <cstring>  // memset


auto bloomflex_patterns_generate(struct bloomflex_s * bloom_filter) -> void
{
  static constexpr unsigned int max_range {63};  // i & max_range = cap values to 63 max
  for(auto i = 0U; i < bloom_filter->pattern_count; i++)
    {
      uint64_t pattern {0};
      for(auto j = 0U; j < bloom_filter->pattern_k; j++)
        {
          uint64_t onebit = 1ULL << (rand_64() & max_range);  // 0 <= shift <= 63
          while ((pattern & onebit) != 0U) {
            onebit = 1ULL << (rand_64() & max_range);
          }
          pattern |= onebit;
        }
      bloom_filter->patterns[i] = pattern;
    }
}


auto bloomflex_init(const uint64_t size, const unsigned int n_hash_functions) -> struct bloomflex_s *
{
  /* Input size is in bytes for full bitmap */

  static constexpr unsigned int multiplier {16};  // multiply by 65,536
  static constexpr unsigned int divider {3};  // divide by 8

  auto * bloom_filter = static_cast<struct bloomflex_s *>(xmalloc(sizeof(struct bloomflex_s)));
  bloom_filter->size = size >> divider;  // divide by 8 to get number of uint64

  bloom_filter->pattern_shift = multiplier;
  bloom_filter->pattern_count = 1U << bloom_filter->pattern_shift;
  bloom_filter->pattern_mask = bloom_filter->pattern_count - 1;
  bloom_filter->pattern_k = n_hash_functions;

  bloom_filter->patterns = static_cast<uint64_t *>(xmalloc(bloom_filter->pattern_count * sizeof(uint64_t)));
  bloomflex_patterns_generate(bloom_filter);

  // b->bitmap = new uint64_t[size / sizeof(uint64_t)];  fix for std::bad_alloc crash
  bloom_filter->bitmap = static_cast<uint64_t *>(xmalloc(size));
  std::memset(bloom_filter->bitmap, UINT8_MAX, size);

  return bloom_filter;
}


auto bloomflex_exit(struct bloomflex_s * bloom_filter) -> void
{
  xfree(bloom_filter->bitmap);
  xfree(bloom_filter->patterns);
  xfree(bloom_filter);
}
