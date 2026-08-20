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

#include "search_dispatch.hpp"
#include "ssse3.hpp"  // dprofile_shuffle8/16 (only declared when __SSSE3__)
#include "sse41.hpp"  // align_cells_*_16_sse41 (only declared when __SSE4_1__)
#include <emmintrin.h>  // _mm_set_epi8/16 (SSE2)
#include "../../utils/cpu_features.hpp"  // Cpu_features
#include "../../utils/search_data.hpp"  // Dprofile_8/16, Score_matrix_8/16, Dseq_8/16
#include <cstdint>  // uint64_t


auto make_T0_16() -> VECTORTYPE16 {
  return _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, -1);
}

auto make_T0_8() -> VECTORTYPE8 {
  return _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0,
                      0, 0, 0, 0, 0, 0, 0, -1);
}


auto dispatch_dprofile16(Cpu_features const & cpu_features,
                         Dprofile_16 & dprofile_a,
                         Score_matrix_16 const & score_matrix_a,
                         Dseq_16 const & dseq_a) -> void
{
#ifdef __SSSE3__
  if (cpu_features.ssse3) {
    dprofile_shuffle16(dprofile_a, score_matrix_a, dseq_a);
    return;
  }
#else
  static_cast<void>(cpu_features);
#endif
  dprofile_fill16(dprofile_a, score_matrix_a, dseq_a);
}

auto dispatch_dprofile8(Cpu_features const & cpu_features,
                        Dprofile_8 & dprofile_a,
                        Score_matrix_8 const & score_matrix_a,
                        Dseq_8 const & dseq_a) -> void
{
#ifdef __SSSE3__
  if (cpu_features.ssse3) {
    dprofile_shuffle8(dprofile_a, score_matrix_a, dseq_a);
    return;
  }
#else
  static_cast<void>(cpu_features);
#endif
  dprofile_fill8(dprofile_a, score_matrix_a, dseq_a);
}


auto dispatch_align_regular_16(Cpu_features const & cpu_features,
                               VECTORTYPE16 * const Sm, VECTORTYPE16 * const hep,
                               VECTORTYPE16 ** const qp,
                               VECTORTYPE16 const & Qm, VECTORTYPE16 const & Rm,
                               uint64_t const ql, VECTORTYPE16 const & F0,
                               uint64_t * const dir_long, VECTORTYPE16 const & H0) -> void
{
#ifdef __SSE4_1__
  if (cpu_features.sse41) {
    align_cells_regular_16_sse41(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0);
    return;
  }
#else
  static_cast<void>(cpu_features);
#endif
  align_cells_regular_16(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0);
}

auto dispatch_align_masked_16(Cpu_features const & cpu_features,
                              VECTORTYPE16 * const Sm, VECTORTYPE16 * const hep,
                              VECTORTYPE16 ** const qp,
                              VECTORTYPE16 const & Qm, VECTORTYPE16 const & Rm,
                              uint64_t const ql, VECTORTYPE16 const & F0,
                              uint64_t * const dir_long, VECTORTYPE16 const & H0,
                              VECTORTYPE16 const * const Mm, VECTORTYPE16 * const MQ,
                              VECTORTYPE16 const * const MR, VECTORTYPE16 const * const MQ0) -> void
{
#ifdef __SSE4_1__
  if (cpu_features.sse41) {
    align_cells_masked_16_sse41(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, Mm, MQ, MR, MQ0);
    return;
  }
#else
  static_cast<void>(cpu_features);
#endif
  align_cells_masked_16(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, Mm, MQ, MR, MQ0);
}
