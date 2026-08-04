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

#ifndef SWARM_ARCH_X86_64_SSE41_H
#define SWARM_ARCH_X86_64_SSE41_H


#ifdef __SSE4_1__

/*
  SSE4.1 specific code for x86-64

  Only include if __SSE4_1__ is defined, which is done by the
  gcc compiler when the -msse4.1 option or similar is given.

  This code requires the _mm_min_epu16 intrinsic implemented
  with the PMINUW instruction on the CPU. That instruction was
  available starting with the Penryn architecture in 2008.
*/

#include <cstdint>
#include <smmintrin.h>


using VECTORTYPE = __m128i;

auto align_cells_regular_16_sse41(VECTORTYPE * Sm,
                                  VECTORTYPE * hep,
                                  VECTORTYPE ** qp,
                                  VECTORTYPE const * Qm,
                                  VECTORTYPE const * Rm,
                                  uint64_t ql,
                                  VECTORTYPE const * F0,
                                  uint64_t * dir_long,
                                  VECTORTYPE const * H0) -> void;

auto align_cells_masked_16_sse41(VECTORTYPE * Sm,
                                 VECTORTYPE * hep,
                                 VECTORTYPE ** qp,
                                 VECTORTYPE const * Qm,
                                 VECTORTYPE const * Rm,
                                 uint64_t ql,
                                 VECTORTYPE const * F0,
                                 uint64_t * dir_long,
                                 VECTORTYPE const * H0,
                                 VECTORTYPE const * Mm,
                                 VECTORTYPE * MQ,
                                 VECTORTYPE const * MR,
                                 VECTORTYPE const * MQ0) -> void;

#endif

#endif  // SWARM_ARCH_X86_64_SSE41_H
