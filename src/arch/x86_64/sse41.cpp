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
#include <cstddef>  // std::ptrdiff_t
#include <limits>  // std::numeric_limits
// C++17 refactoring: [[maybe_unused]]
constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
#endif

using WORD = uint16_t;

namespace {

// The lane operations and direction-word type for this file's width, handed
// to the shared kernel (utils/align_cells.hpp). min() is one PMINUW, where the baseline
// v_min16 has to emulate it; that is the whole of what distinguishes this
// translation unit's kernel from search16.cpp's.
struct Ops_sse41 {
  using Dir_word = WORD;
  static auto add(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_add16(lhs, rhs); }
  static auto sub(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_sub16(lhs, rhs); }
  static auto min(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return _mm_min_epu16(lhs, rhs); }
  static auto mask_eq(VECTORTYPE const lhs, VECTORTYPE const rhs) -> Dir_word { return v_mask_eq16(lhs, rhs); }
  static auto zero() -> VECTORTYPE { return v_zero16(); }
};






// Last, inside this anonymous namespace: see the note at the same point in
// search16.cpp, and utils/align_cells.hpp.
#include "../../utils/align_cells.hpp"



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
  align_cells<Ops_sse41>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, no_mask);
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
  align_cells<Ops_sse41>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, masks);
  *MQ = masks.mq;
}

#else
#error __SSE4_1__ not defined
#endif
