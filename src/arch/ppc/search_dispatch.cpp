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
#include <altivec.h>
#include <cstdint>  // uint64_t
#include <limits>


auto make_T0_16() -> VECTORTYPE16 {
  static constexpr auto unsigned_short_max = std::numeric_limits<unsigned short>::max();
  const VECTORTYPE16 result = { unsigned_short_max, 0, 0, 0, 0, 0, 0, 0 };
  return result;
}

auto make_T0_8() -> VECTORTYPE8 {
  static constexpr auto uchar_max = std::numeric_limits<unsigned char>::max();
  const VECTORTYPE8 result = { uchar_max, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0 };
  return result;
}


// ppc64le has a single Altivec/VSX code path: forward to the generic
// kernels and ignore cpu_features.
auto dispatch_dprofile16(Cpu_features const & cpu_features,
                         Dprofile_16 & dprofile_a,
                         Score_matrix_16 const & score_matrix_a,
                         Dseq_16 const & dseq_a) -> void
{
  static_cast<void>(cpu_features);
  dprofile_fill16(dprofile_a, score_matrix_a, dseq_a);
}

auto dispatch_dprofile8(Cpu_features const & cpu_features,
                        Dprofile_8 & dprofile_a,
                        Score_matrix_8 const & score_matrix_a,
                        Dseq_8 const & dseq_a) -> void
{
  static_cast<void>(cpu_features);
  dprofile_fill8(dprofile_a, score_matrix_a, dseq_a);
}


auto dispatch_align_regular_16(Cpu_features const & cpu_features,
                               VECTORTYPE16 * const Sm, VECTORTYPE16 * const hep,
                               VECTORTYPE16 ** const qp,
                               VECTORTYPE16 const & Qm, VECTORTYPE16 const & Rm,
                               uint64_t const ql, VECTORTYPE16 const & F0,
                               uint64_t * const dir_long, VECTORTYPE16 const & H0) -> void
{
  static_cast<void>(cpu_features);
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
  static_cast<void>(cpu_features);
  align_cells_masked_16(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, Mm, MQ, MR, MQ0);
}
