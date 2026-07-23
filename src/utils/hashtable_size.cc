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

#include "hashtable_size.h"
#include <algorithm>  // std::max
#include <cassert>
#include <cstdint>
#include <cmath>


auto compute_hashtable_size(const uint64_t sequence_count) -> uint64_t {
  // adjust hash table size for at most 70% fill rate (7/10th);
  // i.e. calculate the smallest power of two not smaller than
  // 10/7 times the number of sequences.
  // Note that hash table size can be at least 2^1 and at most 2^63.
  // C++20: refactor with std::bit_ceil()
  static constexpr unsigned int numerator {7};
  static constexpr unsigned int denominator {10};
  static constexpr double divider {2.0};
  static_assert(numerator != 0, "Error: will result in a divide-by-zero");
  assert(sequence_count < 6456360425798343065); // (7 * 2^63 / 10) otherwise hashtable_size > 2^63
  // Scale in floating point: the same product computed in uint64_t,
  // denominator * (sequence_count + 1), overflows for sequence_count
  // above ~1.8e18, i.e. below the assert bound above.
  auto const scaled =
    static_cast<double>(denominator) * (static_cast<double>(sequence_count) + 1.0) / numerator;
  // GCC 6 to 9: std::log2 is not a member of std! (replace with log(x) / log(2.0) for now)
  auto const size = static_cast<uint64_t>(std::pow(divider, std::ceil(std::log(scaled) / std::log(divider))));
  return std::max<uint64_t>(2, size);  // at least 2^1, as documented above
}


/* old function

auto compute_hashtable_size(const uint64_t sequence_count) -> uint64_t {
  static constexpr unsigned int numerator {7};
  static constexpr unsigned int denominator {10};
  uint64_t hashtablesize = 2;
  while (denominator * sequence_count > numerator * hashtablesize) {
    hashtablesize *= 2;
  }
  return hashtablesize;
}

Note that the new function yields hashtablesize values twice bigger
than the values produced by the old function for the following number
of sequences (n), for 1 < n < 10^9:

n       old     new
------------------
11	16	32
179	256	512
2867	4096	8192
45875	65536	131072
734003	1048576	2097152
11744051	16777216	33554432
187904819	268435456	536870912

For all other values (appr. a billion), both functions yield the same
results.

Differences are more frequent for low n values. Adding a small amount
to the sequence count '(sequence_count + 1)' helps reduce the number
of discrepencies to the list above.

I don't think the new function will change swarm's results or
performances.

*/
