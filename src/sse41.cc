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


#ifdef __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#include "utils/intrinsics_to_functions_x86_64.h"
#endif

#ifdef __SSE4_1__

/*
  SSE4.1 specific code for x86-64

  Only include if __SSE4_1__ is defined, which is done by the
  gcc compiler when the -msse4.1 option or similar is given.

  This code requires the _mm_min_epu16 intrinsic implemented
  with the PMINUW instruction on the CPU. That instruction was
  available starting with the Penryn architecture in 2008.
*/

#include <cstdint>  //uint64_t
#include <smmintrin.h>  // _mm_min_epu16

using VECTORTYPE = __m128i;
using WORD = unsigned short;

// refactoring: v_min16 exists and is more complicated
auto v_min(VECTORTYPE lhs, VECTORTYPE rhs) -> VECTORTYPE {
  return _mm_min_epu16(lhs, rhs);
}


inline void onestep_16_sse41(VECTORTYPE & H,
                             VECTORTYPE & N,
                             VECTORTYPE & F,
                             VECTORTYPE V,
                             unsigned short * DIR,
                             VECTORTYPE & E,
                             VECTORTYPE QR,
                             VECTORTYPE R)
{
  H = v_add16(H, V);
  const auto W = H;
  H = v_min(H, F);
  *(DIR + 0) = v_mask_eq16(W, H);
  H = v_min(H, E);
  *(DIR + 1) = v_mask_eq16(H, E);
  N = H;
  H = v_add16(H, QR);
  F = v_add16(F, R);
  E = v_add16(E, R);
  F = v_min(H, F);
  *(DIR + 2) = v_mask_eq16(H, F);
  E = v_min(H, E);
  *(DIR + 3) = v_mask_eq16(H, E);
}


void align_cells_regular_16_sse41(VECTORTYPE * Sm,
                                  VECTORTYPE * hep,
                                  VECTORTYPE ** qp,
                                  VECTORTYPE * Qm,
                                  VECTORTYPE * Rm,
                                  uint64_t ql,
                                  VECTORTYPE * F0,
                                  uint64_t * dir_long,
                                  VECTORTYPE * H0)
{
  static constexpr auto step = 16U;
  static constexpr auto offset0 = 0U;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  const auto Q = *Qm;
  const auto R = *Rm;

  auto f0 = *F0;
  auto f1 = v_add16(f0, R);
  auto f2 = v_add16(f1, R);
  auto f3 = v_add16(f2, R);

  auto h0 = *H0;
  auto h1 = v_sub16(f0, Q);
  auto h2 = v_add16(h1, R);
  auto h3 = v_add16(h2, R);

  auto h5 = v_zero16();
  auto h6 = v_zero16();
  auto h7 = v_zero16();
  auto h8 = v_zero16();

  for(auto i = 0ULL; i < ql; i++)
    {
      VECTORTYPE *x = qp[i + 0];
      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];
      onestep_16_sse41(h0, h5, f0, x[0], dir + step * i + offset0, E, Q, R);
      onestep_16_sse41(h1, h6, f1, x[1], dir + step * i + offset1, E, Q, R);
      onestep_16_sse41(h2, h7, f2, x[2], dir + step * i + offset2, E, Q, R);
      onestep_16_sse41(h3, h8, f3, x[3], dir + step * i + offset3, E, Q, R);
      hep[2 * i + 0] = h8;
      hep[2 * i + 1] = E;
      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}


void align_cells_masked_16_sse41(VECTORTYPE * Sm,
                                 VECTORTYPE * hep,
                                 VECTORTYPE ** qp,
                                 VECTORTYPE * Qm,
                                 VECTORTYPE * Rm,
                                 uint64_t ql,
                                 VECTORTYPE * F0,
                                 uint64_t * dir_long,
                                 VECTORTYPE * H0,
                                 VECTORTYPE * Mm,
                                 VECTORTYPE * MQ,
                                 VECTORTYPE * MR,
                                 VECTORTYPE * MQ0)
{
  static constexpr auto step = 16U;
  static constexpr auto offset0 = 0U;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  const auto Q = *Qm;
  const auto R = *Rm;

  auto f0 = *F0;
  auto f1 = v_add16(f0, R);
  auto f2 = v_add16(f1, R);
  auto f3 = v_add16(f2, R);

  auto h0 = *H0;
  auto h1 = v_sub16(f0, Q);
  auto h2 = v_add16(h1, R);
  auto h3 = v_add16(h2, R);

  auto h5 = v_zero16();
  auto h6 = v_zero16();
  auto h7 = v_zero16();
  auto h8 = v_zero16();

  for(auto i = 0ULL; i < ql; i++)
    {
      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];
      VECTORTYPE *x = qp[i + 0];

      /* mask h4 and E */
      h4 = v_sub16(h4, *Mm);
      E  = v_sub16(E,  *Mm);

      /* init h4 and E */
      h4 = v_add16(h4, *MQ);
      E  = v_add16(E,  *MQ);
      E  = v_add16(E,  *MQ0);

      /* update MQ */
      *MQ = v_add16(*MQ,  *MR);


      onestep_16_sse41(h0, h5, f0, x[0], dir + step * i + offset0, E, Q, R);
      onestep_16_sse41(h1, h6, f1, x[1], dir + step * i + offset1, E, Q, R);
      onestep_16_sse41(h2, h7, f2, x[2], dir + step * i + offset2, E, Q, R);
      onestep_16_sse41(h3, h8, f3, x[3], dir + step * i + offset3, E, Q, R);
      hep[2 * i + 0] = h8;
      hep[2 * i + 1] = E;

      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}

#else
#error __SSE4_1__ not defined
#endif
#endif
