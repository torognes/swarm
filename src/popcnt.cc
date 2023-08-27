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


#ifdef __x86_64__
#ifdef __POPCNT__

#include "swarm.h"
#include <popcntintrin.h>
#include <cstdint>  // uint64_t


/*
  POPCNT specific code for x86-64

  Only include if __POPCNT__ is defined, which is done by the
  gcc compiler when the -mpopcnt option or similar is given.

  This code requires the _mm_popcnt_u64 intrinsic implemented
  with the POPCNT instruction on the CPU. That instruction was
  available starting with the Nahalem architecture in 2008.
*/

auto compareqgramvectors_popcnt(unsigned char * qgram_a, unsigned char * qgram_b) -> uint64_t
{
  /* Count number of different bits */
  /* Uses the POPCNT instruction, requires CPU with this feature */

  auto *ap = reinterpret_cast<uint64_t*>(qgram_a);
  auto *bp = reinterpret_cast<uint64_t*>(qgram_b);
  uint64_t count {0};

  while (reinterpret_cast<unsigned char*>(ap) < qgram_a + qgramvectorbytes) {
    count += _mm_popcnt_u64(*ap++ ^ *bp++);
  }

  return count;
}

#else
#error __POPCNT__ not defined
#endif
#endif
