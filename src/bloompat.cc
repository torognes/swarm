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
#include <algorithm>  // std::max(), std::fill()
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // uint64_t, uint8_t
#include <cstring>
#include <iterator>  // std::next
#include <limits>


// --------------------------------------------------- used only in bloompat.cc

auto bloom_adr(struct bloom_s * bloom_filter, const uint64_t hash) -> uint64_t *
{
  auto const position = static_cast<std::ptrdiff_t>((hash >> bloom_pattern_shift) & bloom_filter->mask);
  return std::next(bloom_filter->bitmap, position);
}

auto bloom_pat(struct bloom_s * bloom_filter, uint64_t hash) -> uint64_t
{
  return bloom_filter->patterns[hash & bloom_pattern_mask];
}

// ------------------------------------------------------------------------ end

// used in algod1.cc
auto bloom_set(struct bloom_s * bloom_filter, uint64_t hash) -> void
{
  *bloom_adr(bloom_filter, hash) &= compl bloom_pat(bloom_filter, hash);
}

// used in algod1.cc
auto bloom_get(struct bloom_s * bloom_filter, uint64_t hash) -> bool
{
  return (*bloom_adr(bloom_filter, hash) & bloom_pat(bloom_filter, hash)) == 0U;
}


auto bloom_patterns_generate(struct bloom_s & bloom_filter) -> void
{
  static constexpr auto max_range = 63U;  // i & max_range = cap values to 63 max
  static constexpr auto n_loops = 8U;
  for(auto & pattern : bloom_filter.patterns)
    {
      pattern = 0;
      for(auto i = 0U; i < n_loops; ++i)
        {
          uint64_t onebit = 1ULL << (rand_64() & max_range);  // 0 <= shift <= 63
          while ((pattern & onebit) != 0) {
            onebit = 1ULL << (rand_64() & max_range);
          }
          pattern |= onebit;
        }
    }
}


auto bloom_zap(struct bloom_s & bloom_filter) -> void
{
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
  std::fill(bloom_filter.bitmap_v.begin(), bloom_filter.bitmap_v.end(), uint8_max);
}


auto bloom_init(uint64_t size, struct bloom_s & bloom_filter) -> struct bloom_s *
{
  // Size is in bytes for full bitmap, must be power of 2,
  // at least 8
  // assert(std::has_single_bit(size));  // C++20 refactoring: is power of 2?
  static constexpr uint64_t bytes_per_uint64 {8};
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();

  size = std::max(size, bytes_per_uint64);

  bloom_filter.size = size;

  bloom_filter.mask = (size >> 3U) - 1;

  bloom_filter.bitmap_v.resize(size, uint8_max);
  bloom_filter.bitmap = bloom_filter.bitmap_v.data();

  bloom_patterns_generate(bloom_filter);

  return &bloom_filter;
}


auto bloom_exit(struct bloom_s * bloom_filter) -> void
{
  bloom_filter->bitmap = nullptr;
  bloom_filter = nullptr;
}
