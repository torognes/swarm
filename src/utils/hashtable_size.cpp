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

#include "hashtable_size.hpp"
#include <cassert>
#include <cstdint>


auto compute_hashtable_size(const uint64_t sequence_count) -> uint64_t {
  // adjust hash table size for at most 70% fill rate (7/10th);
  // i.e. calculate the smallest power of two not smaller than
  // 10/7 times the number of sequences.
  // Note that hash table size can be at least 2^1 and at most 2^63.
  // C++20: refactor with std::bit_ceil()
  static constexpr uint64_t numerator {7};
  static constexpr uint64_t denominator {10};
  static constexpr uint64_t smallest {2};                   // 2^1, as documented above
  static constexpr uint64_t largest {uint64_t{1} << 63};    // 2^63, likewise
  static_assert(numerator != 0, "Error: will result in a divide-by-zero");
  assert(sequence_count < 6456360425798343065); // (7 * 2^63 / 10) otherwise hashtable_size > 2^63

  // Integer arithmetic, exactly. This used to scale in double and take
  // pow(2, ceil(log(scaled) / log(2))), which rounds four times; computing
  // log2 as log(x)/log(2) is the one that costs correctness, because near an
  // exact power of two the quotient lands on the wrong side and ceil then
  // picks the wrong exponent. It errs in both directions: at
  // sequence_count = 394064967394918 the old form returned 2^49 where the
  // answer is 2^50, and at 1576259869579672 it returned 2^52 where the
  // answer is 2^51.
  //
  // The condition wanted is 'size >= 10 * (sequence_count + 1) / 7' with a
  // real division, i.e. size >= ceil(10 * m / 7). The product 10 * m
  // overflows uint64_t inside the range asserted above, which is why the
  // double was there; it is never formed here. Writing m as 7q + r gives
  // 10m/7 == 10q + 10r/7 with r in [0, 6], so
  //   ceil(10m/7) == 10q + ceil(10r/7) == 10q + (10r + 6) / 7
  // whose largest term is 10q <= 2^63 for any m the assert admits.
  auto const scaled = sequence_count + 1;
  auto const threshold = (denominator * (scaled / numerator))
                       + (((denominator * (scaled % numerator)) + numerator - 1) / numerator);

  // The cap is what keeps this total rather than merely asserted: without
  // it, a threshold above 2^63 would double past it and wrap to zero,
  // looping forever in a build with NDEBUG. Returning 2^63 there is also
  // what the previous implementation did.
  auto size = smallest;
  while ((size < threshold) and (size < largest)) {
    size *= 2;
  }
  return size;
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
