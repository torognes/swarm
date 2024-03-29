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

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
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
#define v_min(a, b) _mm_min_epu16((a), (b))
#define v_add(a, b) _mm_adds_epu16((a), (b))
#define v_sub(a, b) _mm_subs_epu16((a), (b))
#define v_dup(a) _mm_set1_epi16(a)
#define v_zero v_dup(0)
#define v_mask_eq(a, b) static_cast<unsigned short>(_mm_movemask_epi8(_mm_cmpeq_epi16((a), (b))))

inline void onestep_16_sse41(VECTORTYPE & H,
                             VECTORTYPE & N,
                             VECTORTYPE & F,
                             VECTORTYPE V,
                             unsigned short * DIR,
                             VECTORTYPE & E,
                             VECTORTYPE QR,
                             VECTORTYPE R)
{
  VECTORTYPE W;

  H = v_add(H, V);
  W = H;
  H = v_min(H, F);
  *(DIR+0) = v_mask_eq(W, H);
  H = v_min(H, E);
  *(DIR+1) = v_mask_eq(H, E);
  N = H;
  H = v_add(H, QR);
  F = v_add(F, R);
  E = v_add(E, R);
  F = v_min(H, F);
  *(DIR+2) = v_mask_eq(H, F);
  E = v_min(H, E);
  *(DIR+3) = v_mask_eq(H, E);
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
  static constexpr unsigned int step {16};
  static constexpr unsigned int offset0 {0};
  static constexpr unsigned int offset1 {offset0 + 4};
  static constexpr unsigned int offset2 {offset1 + 4};
  static constexpr unsigned int offset3 {offset2 + 4};

  VECTORTYPE Q;
  VECTORTYPE R;
  VECTORTYPE E;
  VECTORTYPE h0;
  VECTORTYPE h1;
  VECTORTYPE h2;
  VECTORTYPE h3;
  VECTORTYPE h4;
  VECTORTYPE h5;
  VECTORTYPE h6;
  VECTORTYPE h7;
  VECTORTYPE h8;
  VECTORTYPE f0;
  VECTORTYPE f1;
  VECTORTYPE f2;
  VECTORTYPE f3;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  Q = *Qm;
  R = *Rm;

  f0 = *F0;
  f1 = v_add(f0, R);
  f2 = v_add(f1, R);
  f3 = v_add(f2, R);

  h0 = *H0;
  h1 = v_sub(f0, Q);
  h2 = v_add(h1, R);
  h3 = v_add(h2, R);

  h5 = v_zero;
  h6 = v_zero;
  h7 = v_zero;
  h8 = v_zero;

  for(auto i = 0ULL; i < ql; i++)
    {
      VECTORTYPE *x {nullptr};

      x = qp[i + 0];
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
  static constexpr unsigned int step {16};
  static constexpr unsigned int offset0 {0};
  static constexpr unsigned int offset1 {offset0 + 4};
  static constexpr unsigned int offset2 {offset1 + 4};
  static constexpr unsigned int offset3 {offset2 + 4};

  VECTORTYPE Q;
  VECTORTYPE R;
  VECTORTYPE E;
  VECTORTYPE h0;
  VECTORTYPE h1;
  VECTORTYPE h2;
  VECTORTYPE h3;
  VECTORTYPE h4;
  VECTORTYPE h5;
  VECTORTYPE h6;
  VECTORTYPE h7;
  VECTORTYPE h8;
  VECTORTYPE f0;
  VECTORTYPE f1;
  VECTORTYPE f2;
  VECTORTYPE f3;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  Q = *Qm;
  R = *Rm;

  f0 = *F0;
  f1 = v_add(f0, R);
  f2 = v_add(f1, R);
  f3 = v_add(f2, R);

  h0 = *H0;
  h1 = v_sub(f0, Q);
  h2 = v_add(h1, R);
  h3 = v_add(h2, R);

  h5 = v_zero;
  h6 = v_zero;
  h7 = v_zero;
  h8 = v_zero;

  for(auto i = 0ULL; i < ql; i++)
    {
      VECTORTYPE *x {nullptr};

      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];
      x = qp[i + 0];

      /* mask h4 and E */
      h4 = v_sub(h4, *Mm);
      E  = v_sub(E,  *Mm);

      /* init h4 and E */
      h4 = v_add(h4, *MQ);
      E  = v_add(E,  *MQ);
      E  = v_add(E,  *MQ0);

      /* update MQ */
      *MQ = v_add(*MQ,  *MR);


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
