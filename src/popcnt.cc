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

#include <popcntintrin.h>
#include "utils/qgram_array.h"
#include <cstdint>  // uint64_t


/*
  POPCNT specific code for x86-64

  Only include if __POPCNT__ is defined, which is done by the
  gcc compiler when the -mpopcnt option or similar is given.

  This code requires the _mm_popcnt_u64 intrinsic implemented
  with the POPCNT instruction on the CPU. That instruction was
  available starting with the Nehalem architecture in 2008.
*/

auto compareqgramvectors_popcnt(unsigned char * qgram_a, unsigned char * qgram_b) -> uint64_t
{
  /* Count number of different bits */
  /* Uses the POPCNT instruction, requires CPU with this feature */

  auto *qgram_vector_a_ptr = reinterpret_cast<uint64_t*>(qgram_a);
  auto *qgram_vector_b_ptr = reinterpret_cast<uint64_t*>(qgram_b);
  uint64_t count {0};

  while (reinterpret_cast<unsigned char*>(qgram_vector_a_ptr) < qgram_a + qgramvectorbytes) {
    count += static_cast<uint64_t>(_mm_popcnt_u64(*qgram_vector_a_ptr++ ^ *qgram_vector_b_ptr++));
  }

  return count;
}

#else
#error __POPCNT__ not defined
#endif
#endif
