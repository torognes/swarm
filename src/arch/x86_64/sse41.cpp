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

// The SSE4.1 unsigned 16-bit minimum: one PMINUW, where the baseline
// v_min16 has to emulate it. This is the whole of what distinguishes this
// translation unit's kernel from search16.cpp's -- see
// utils/align_cells_16.hpp.
struct Min_sse41 {
  static auto min(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE {
    return _mm_min_epu16(lhs, rhs);
  }
};




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


// Last, inside this anonymous namespace: see the note at the same point in
// search16.cpp, and utils/align_cells_16.hpp.
#include "../../utils/align_cells_16.hpp"



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
  align_cells_16<Min_sse41>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, no_mask);
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
  align_cells_16<Min_sse41>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, masks);
  *MQ = masks.mq;
}

#else
#error __SSE4_1__ not defined
#endif
