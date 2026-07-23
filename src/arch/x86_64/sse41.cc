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


#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#include "intrinsics_to_functions.h"
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
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
#include "sse41.h"

#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
#include <limits>
constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
#endif

using VECTORTYPE = __m128i;
using WORD = unsigned short;

namespace {

// refactoring: v_min16 exists and is more complicated
auto v_min(VECTORTYPE lhs, VECTORTYPE rhs) -> VECTORTYPE {
  return _mm_min_epu16(lhs, rhs);
}


inline auto onestep_16_sse41(VECTORTYPE & H,
                             VECTORTYPE & N,
                             VECTORTYPE & F,
                             VECTORTYPE V,
                             WORD * DIR,
                             VECTORTYPE & E,
                             VECTORTYPE QR,
                             VECTORTYPE R) -> void
{
  H = v_add16(H, V);
  const auto W = H;
  H = v_min(H, F);
  DIR[0] = v_mask_eq16(W, H);  // subscript, not std::next: hot loop, see align_cells
  H = v_min(H, E);
  DIR[1] = v_mask_eq16(H, E);
  N = H;
  H = v_add16(H, QR);
  F = v_add16(F, R);
  E = v_add16(E, R);
  F = v_min(H, F);
  DIR[2] = v_mask_eq16(H, F);
  E = v_min(H, E);
  DIR[3] = v_mask_eq16(H, E);
}


// One block of cells, shared by the regular and masked SSE4.1 kernels.
// The masked variant differs only by a per-iteration adjustment of h4
// and E (consuming the Mm / MQ / MR / MQ0 vectors); 'masked' is a
// compile-time flag, so the regular instantiation drops that block
// entirely and never dereferences the (null) masking pointers.
template <bool masked>
auto align_cells_16_sse41(VECTORTYPE * Sm,
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
                          VECTORTYPE const * MQ0) -> void
{
  static constexpr auto step = 16;
  static constexpr auto offset0 = 0;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * const dir = reinterpret_cast<WORD *>(dir_long);

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

  assert(ql <= max_ptrdiff);
  assert(ql <= ((max_ptrdiff - 1) / 2));  // max 'E' offset
  assert(ql <= ((max_ptrdiff - offset3) / step));  // max 'dir' offset
  auto const ql_signed = static_cast<std::ptrdiff_t>(ql);
  // Performance: index this hot loop with the subscript operator (and
  // &dir[...] for the 'dir' pointer arguments) rather than std::next().
  // The std::next() form (commit 8c6925f, "pro-bounds-pointer-arithmetic")
  // pessimized this kernel by ~20% on d>1 18SV9. Subscripting restores it and
  // stays clang-tidy clean: pos is signed (no -Wsign-conversion) and operator[]
  // is not flagged by cppcoreguidelines-pro-bounds-pointer-arithmetic.
  for (auto pos = 0LL; pos < ql_signed; ++pos) {
      VECTORTYPE const * x = qp[pos];
      h4 = hep[(2 * pos) + 0];
      E  = hep[(2 * pos) + 1];

      if (masked)
        {
          /* mask h4 and E */
          h4 = v_sub16(h4, *Mm);
          E  = v_sub16(E,  *Mm);

          /* init h4 and E */
          h4 = v_add16(h4, *MQ);
          E  = v_add16(E,  *MQ);
          E  = v_add16(E,  *MQ0);

          /* update MQ */
          *MQ = v_add16(*MQ,  *MR);
        }

      onestep_16_sse41(h0, h5, f0, x[0], &dir[(step * pos) + offset0], E, Q, R);
      onestep_16_sse41(h1, h6, f1, x[1], &dir[(step * pos) + offset1], E, Q, R);
      onestep_16_sse41(h2, h7, f2, x[2], &dir[(step * pos) + offset2], E, Q, R);
      onestep_16_sse41(h3, h8, f3, x[3], &dir[(step * pos) + offset3], E, Q, R);
      hep[(2 * pos) + 0] = h8;
      hep[(2 * pos) + 1] = E;
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

}  // namespace


auto align_cells_regular_16_sse41(VECTORTYPE * Sm,
                                  VECTORTYPE * hep,
                                  VECTORTYPE ** qp,
                                  VECTORTYPE const * Qm,
                                  VECTORTYPE const * Rm,
                                  uint64_t ql,
                                  VECTORTYPE const * F0,
                                  uint64_t * dir_long,
                                  VECTORTYPE const * H0) -> void
{
  align_cells_16_sse41<false>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0,
                              nullptr, nullptr, nullptr, nullptr);
}


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
                                 VECTORTYPE const * MQ0) -> void
{
  align_cells_16_sse41<true>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, Mm, MQ, MR, MQ0);
}

#else
#error __SSE4_1__ not defined
#endif
