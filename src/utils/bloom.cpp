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

#include "bloom.hpp"
#include "pseudo_rng.hpp"
#include <cassert>
#include <cstdint>  // uint64_t
#include <vector>


namespace bloom_detail {

auto generate_patterns(std::vector<uint64_t> & patterns,
                       uint64_t const pattern_k) -> void {
  static constexpr auto max_range = 63U;  // i & max_range = cap values to 63 max
  for (auto & pattern : patterns) {
    assert(pattern == 0);  // value-initialized by the vector constructor
    for (auto j = 0U; j < pattern_k; ++j) {
      uint64_t onebit = 1ULL << (rand_64() & max_range);  // 0 <= shift <= 63
      while ((pattern & onebit) != 0U) {
        onebit = 1ULL << (rand_64() & max_range);
      }
      pattern |= onebit;
    }
  }
}

}  // namespace bloom_detail
