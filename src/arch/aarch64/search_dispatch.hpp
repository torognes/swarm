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

#ifndef SWARM_ARCH_AARCH64_SEARCH_DISPATCH_H
#define SWARM_ARCH_AARCH64_SEARCH_DISPATCH_H

#include "../../utils/cpu_features.hpp"  // Cpu_features
#include "../../utils/search_data.hpp"  // BYTE, WORD, Score_matrix_8/16, Dseq_8/16
#include <arm_neon.h>  // uint16x8_t, uint8x16_t
#include <cstdint>  // uint64_t
#include <vector>  // std::vector


// Per-architecture entry points for search8 / search16. AArch64 has a
// single NEON code path, so these only build the diagonal mask seed (T0)
// and forward to the generic kernels; cpu_features is ignored. The
// interface matches the x86_64 dispatcher so the search functions stay
// architecture-agnostic. Mirrors the scheme used for compareqgramvectors
// (see utils/qgram_compare.hpp).

using VECTORTYPE16 = uint16x8_t;
using VECTORTYPE8 = uint8x16_t;


// Initial diagonal mask: highest channel set, all others zero.
auto make_T0_16() -> VECTORTYPE16;
auto make_T0_8() -> VECTORTYPE8;


auto dispatch_dprofile16(Cpu_features const & cpu_features,
                         std::vector<WORD> & dprofile_v,
                         Score_matrix_16 const & score_matrix_a,
                         Dseq_16 const & dseq_a) -> void;
auto dispatch_dprofile8(Cpu_features const & cpu_features,
                        std::vector<BYTE> & dprofile_v,
                        Score_matrix_8 const & score_matrix_a,
                        Dseq_8 const & dseq_a) -> void;


auto dispatch_align_regular_16(Cpu_features const & cpu_features,
                               VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                               VECTORTYPE16 ** qp,
                               VECTORTYPE16 const & Qm, VECTORTYPE16 const & Rm,
                               uint64_t ql, VECTORTYPE16 const & F0,
                               uint64_t * dir_long, VECTORTYPE16 const & H0) -> void;
auto dispatch_align_masked_16(Cpu_features const & cpu_features,
                              VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                              VECTORTYPE16 ** qp,
                              VECTORTYPE16 const & Qm, VECTORTYPE16 const & Rm,
                              uint64_t ql, VECTORTYPE16 const & F0,
                              uint64_t * dir_long, VECTORTYPE16 const & H0,
                              VECTORTYPE16 const * Mm, VECTORTYPE16 * MQ,
                              VECTORTYPE16 const * MR, VECTORTYPE16 const * MQ0) -> void;


// Generic NEON kernels, defined in search16.cpp / search8.cpp and called by
// the dispatchers above.
auto dprofile_fill16(std::vector<WORD> & dprofile_v,
                     Score_matrix_16 const & score_matrix_a,
                     Dseq_16 const & dseq_a) -> void;
auto dprofile_fill8(std::vector<BYTE> & dprofile_v,
                    Score_matrix_8 const & score_matrix_a,
                    Dseq_8 const & dseq_a) -> void;
auto align_cells_regular_16(VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                            VECTORTYPE16 ** qp,
                            VECTORTYPE16 const & Qm, VECTORTYPE16 const & Rm,
                            uint64_t ql, VECTORTYPE16 const & F0,
                            uint64_t * dir_long, VECTORTYPE16 const & H0) -> void;
auto align_cells_masked_16(VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                           VECTORTYPE16 ** qp,
                           VECTORTYPE16 const & Qm, VECTORTYPE16 const & Rm,
                           uint64_t ql, VECTORTYPE16 const & F0,
                           uint64_t * dir_long, VECTORTYPE16 const & H0,
                           VECTORTYPE16 const * Mm, VECTORTYPE16 * MQ,
                           VECTORTYPE16 const * MR, VECTORTYPE16 const * MQ0) -> void;

#endif
