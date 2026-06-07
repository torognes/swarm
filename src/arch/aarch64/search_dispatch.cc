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

#include "search_dispatch.h"
#include <arm_neon.h>
#include <cstdint>  // uint16_t, uint8_t
#include <limits>


auto make_T0_16() -> VECTORTYPE16 {
  static constexpr auto uint16_max = std::numeric_limits<uint16_t>::max();
  const VECTORTYPE16 result = { uint16_max, 0, 0, 0, 0, 0, 0, 0 };
  return result;
}

auto make_T0_8() -> VECTORTYPE8 {
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
  const VECTORTYPE8 result = { uint8_max, 0, 0, 0, 0, 0, 0, 0,
                                       0, 0, 0, 0, 0, 0, 0, 0 };
  return result;
}


// AArch64 has a single NEON code path: forward to the generic kernels and
// ignore cpu_features.
auto dispatch_dprofile16(Cpu_features const & cpu_features,
                         unsigned short * dprofile,
                         unsigned short const * score_matrix,
                         unsigned char const * dseq) -> void
{
  static_cast<void>(cpu_features);
  dprofile_fill16(dprofile, score_matrix, dseq);
}

auto dispatch_dprofile8(Cpu_features const & cpu_features,
                        unsigned char * dprofile,
                        unsigned char const * score_matrix,
                        unsigned char const * dseq) -> void
{
  static_cast<void>(cpu_features);
  dprofile_fill8(dprofile, score_matrix, dseq);
}


auto dispatch_align_regular_16(Cpu_features const & cpu_features,
                               VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                               VECTORTYPE16 ** qp,
                               VECTORTYPE16 const * Qm, VECTORTYPE16 const * Rm,
                               uint64_t ql, VECTORTYPE16 const * F0,
                               uint64_t * dir_long, VECTORTYPE16 const * H0) -> void
{
  static_cast<void>(cpu_features);
  align_cells_regular_16(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0);
}

auto dispatch_align_masked_16(Cpu_features const & cpu_features,
                              VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                              VECTORTYPE16 ** qp,
                              VECTORTYPE16 const * Qm, VECTORTYPE16 const * Rm,
                              uint64_t ql, VECTORTYPE16 const * F0,
                              uint64_t * dir_long, VECTORTYPE16 const * H0,
                              VECTORTYPE16 const * Mm, VECTORTYPE16 * MQ,
                              VECTORTYPE16 const * MR, VECTORTYPE16 const * MQ0) -> void
{
  static_cast<void>(cpu_features);
  align_cells_masked_16(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, Mm, MQ, MR, MQ0);
}
