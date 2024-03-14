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

#include "bloomflex.h"
#include "utils/pseudo_rng.h"
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // uint64_t
#include <cstring>  // memset
#include <iterator>  // std::next
#include <limits>


auto bloomflex_adr(struct bloomflex_s * bloom_filter, const uint64_t hash) -> uint64_t *
{
  auto const position = (hash >> bloom_filter->pattern_shift) % bloom_filter->size;
  assert(position <= std::numeric_limits<std::ptrdiff_t>::max());
  auto const signed_position = static_cast<std::ptrdiff_t>(position);
  return std::next(bloom_filter->bitmap, signed_position);
}


auto bloomflex_pat(struct bloomflex_s * bloom_filter, const uint64_t hash) -> uint64_t
{
  auto const position = hash & bloom_filter->pattern_mask;
  assert(position <= std::numeric_limits<std::ptrdiff_t>::max());
  auto const signed_position = static_cast<std::ptrdiff_t>(position);
  return *std::next(bloom_filter->patterns, signed_position);
}


auto bloomflex_set(struct bloomflex_s * bloom_filter, uint64_t hash) -> void
{
  *bloomflex_adr(bloom_filter, hash) &= compl bloomflex_pat(bloom_filter, hash);
}


auto bloomflex_get(struct bloomflex_s * bloom_filter, uint64_t hash) -> bool
{
  return (*bloomflex_adr(bloom_filter, hash) & bloomflex_pat(bloom_filter, hash)) == 0U;
}


auto bloomflex_patterns_generate(struct bloomflex_s & bloom_filter) -> void
{
  static constexpr auto max_range = 63U;  // i & max_range = cap values to 63 max
  for(auto & pattern : bloom_filter.patterns_v)
    {
      pattern = 0;
      for(auto j = 0U; j < bloom_filter.pattern_k; ++j)
        {
          uint64_t onebit = 1ULL << (rand_64() & max_range);  // 0 <= shift <= 63
          while ((pattern & onebit) != 0U) {
            onebit = 1ULL << (rand_64() & max_range);
          }
          pattern |= onebit;
        }
    }
}


auto bloomflex_init(const uint64_t size, const unsigned int n_hash_functions,
                    struct bloomflex_s& bloom_filter) -> struct bloomflex_s *
{
  /* Input size is in bytes for full bitmap */

  static constexpr unsigned int multiplier {16};  // multiply by 65,536
  static constexpr unsigned int divider {3};  // divide by 8
  static constexpr auto uint64_max = std::numeric_limits<uint64_t>::max();

  bloom_filter.size = size >> divider;  // divide by 8 to get number of uint64

  bloom_filter.pattern_shift = multiplier;
  bloom_filter.pattern_count = 1U << bloom_filter.pattern_shift;
  bloom_filter.pattern_mask = bloom_filter.pattern_count - 1;
  bloom_filter.pattern_k = n_hash_functions;

  bloom_filter.patterns_v.resize(bloom_filter.pattern_count);
  bloom_filter.patterns = bloom_filter.patterns_v.data();
  bloomflex_patterns_generate(bloom_filter);

  bloom_filter.bitmap_v.resize(bloom_filter.size, uint64_max);
  bloom_filter.bitmap = bloom_filter.bitmap_v.data();

  return &bloom_filter;
}


auto bloomflex_exit(struct bloomflex_s & bloom_filter) -> void
{
  bloom_filter.bitmap = nullptr;
  bloom_filter.patterns =nullptr;
}
