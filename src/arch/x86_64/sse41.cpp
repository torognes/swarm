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
#include "intrinsics_to_functions.hpp"
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

#include <cstdint>  // uint16_t, uint64_t
#include <smmintrin.h>  // _mm_min_epu16
#include "../../utils/mask_vectors.hpp"  // No_mask, Mask_vectors
#include "sse41.hpp"  // VECTORTYPE

#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
#include <limits>
constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
#endif

using WORD = uint16_t;

namespace {

// refactoring: v_min16 exists and is more complicated
auto v_min(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE {
  return _mm_min_epu16(lhs, rhs);
}


inline auto onestep_16_sse41(VECTORTYPE & H,
                             VECTORTYPE & N,
                             VECTORTYPE & F,
                             VECTORTYPE const V,
                             WORD * const DIR,
                             VECTORTYPE & E,
                             VECTORTYPE const QR,
                             VECTORTYPE const R) -> void
{
  H = v_add16(H, V);
  auto const W = H;
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


// The masking payload: 'mask' selects the channels whose sequence just
// ended, 'mq' is the running gap-open accumulator seeded by the caller,
// 'mr' its per-iteration increment, and 'mq0' the value 'mq' held on
// entry. A plain struct rather than a template on VECTORTYPE: see
// utils/mask_vectors.hpp for why the template form is not usable here.
struct Mask_vectors {
  VECTORTYPE mask;
  VECTORTYPE mq;
  VECTORTYPE mr;
  VECTORTYPE mq0;
};


// The masking step, selected by the type of the kernel's mask argument
// (see utils/mask_vectors.hpp). The No_mask overload is empty, so the
// regular kernel's loop body contains nothing at this point.
inline auto apply_mask(VECTORTYPE & /*h4*/, VECTORTYPE & /*E*/,
                       No_mask const & /*masks*/) -> void
{
}

inline auto apply_mask(VECTORTYPE & h4, VECTORTYPE & E,
                       Mask_vectors & masks) -> void
{
  /* mask h4 and E */
  h4 = v_sub16(h4, masks.mask);
  E  = v_sub16(E,  masks.mask);

  /* init h4 and E */
  h4 = v_add16(h4, masks.mq);
  E  = v_add16(E,  masks.mq);
  E  = v_add16(E,  masks.mq0);

  /* update MQ */
  masks.mq = v_add16(masks.mq,  masks.mr);
}


// One block of cells, shared by the regular and masked SSE4.1 kernels.
// The masked variant differs only by a per-iteration adjustment of h4
// and E; which flavour this is comes from the type of 'masks', so the
// regular instantiation drops that adjustment entirely and is handed no
// masking data at all (see utils/mask_vectors.hpp).
template <typename Masks>
auto align_cells_16_sse41(VECTORTYPE * const Sm,
                          VECTORTYPE * const hep,
                          VECTORTYPE ** const qp,
                          VECTORTYPE const & Qm,
                          VECTORTYPE const & Rm,
                          uint64_t const ql,
                          VECTORTYPE const & F0,
                          uint64_t * const dir_long,
                          VECTORTYPE const & H0,
                          Masks & masks) -> void
{
  static constexpr auto step = 16;
  static constexpr auto offset0 = 0;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * const dir = reinterpret_cast<WORD *>(dir_long);

  auto const Q = Qm;
  auto const R = Rm;

  auto f0 = F0;
  auto f1 = v_add16(f0, R);
  auto f2 = v_add16(f1, R);
  auto f3 = v_add16(f2, R);

  auto h0 = H0;
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
      VECTORTYPE const * const x = qp[pos];
      h4 = hep[(2 * pos) + 0];
      E  = hep[(2 * pos) + 1];

      apply_mask(h4, E, masks);

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


auto align_cells_regular_16_sse41(VECTORTYPE * const Sm,
                                  VECTORTYPE * const hep,
                                  VECTORTYPE ** const qp,
                                  VECTORTYPE const & Qm,
                                  VECTORTYPE const & Rm,
                                  uint64_t const ql,
                                  VECTORTYPE const & F0,
                                  uint64_t * const dir_long,
                                  VECTORTYPE const & H0) -> void
{
  No_mask no_mask;
  align_cells_16_sse41(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, no_mask);
}


auto align_cells_masked_16_sse41(VECTORTYPE * const Sm,
                                 VECTORTYPE * const hep,
                                 VECTORTYPE ** const qp,
                                 VECTORTYPE const & Qm,
                                 VECTORTYPE const & Rm,
                                 uint64_t const ql,
                                 VECTORTYPE const & F0,
                                 uint64_t * const dir_long,
                                 VECTORTYPE const & H0,
                                 VECTORTYPE const * const Mm,
                                 VECTORTYPE * const MQ,
                                 VECTORTYPE const * const MR,
                                 VECTORTYPE const * const MQ0) -> void
{
  Mask_vectors masks {*Mm, *MQ, *MR, *MQ0};
  align_cells_16_sse41(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, masks);
  *MQ = masks.mq;
}

#else
#error __SSE4_1__ not defined
#endif
