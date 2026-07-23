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

#ifndef SWARM_ARCH_X86_64_SEARCH_DISPATCH_H
#define SWARM_ARCH_X86_64_SEARCH_DISPATCH_H

#include "../../utils/cpu_features.hpp"  // Cpu_features
#include <cstdint>  // uint64_t
#include <emmintrin.h>  // __m128i (SSE2)


// Per-architecture entry points for search8 / search16. The diagonal
// mask seed (T0) and the SSSE3 / SSE4.1 run-time dispatch live here, so
// the search functions themselves stay free of architecture and
// instruction-set preprocessor branches.
//
// cpu_features is consulted only when this translation unit is compiled
// with __SSE3__ / __SSE4_1__ defined (see the Makefile per-file flags);
// with the default -march=x86-64 baseline neither is defined and the
// generic SSE2 kernels are always selected. Mirrors the dispatch scheme
// already used for compareqgramvectors (see utils/qgram_compare.hpp).

using VECTORTYPE16 = __m128i;
using VECTORTYPE8 = __m128i;


// Initial diagonal mask: highest channel set, all others zero.
auto make_T0_16() -> VECTORTYPE16;
auto make_T0_8() -> VECTORTYPE8;


// Score-profile construction (SSSE3 shuffle when available, else the
// generic gather).
auto dispatch_dprofile16(Cpu_features const & cpu_features,
                         unsigned short * dprofile,
                         unsigned short const * score_matrix,
                         unsigned char const * dseq) -> void;
auto dispatch_dprofile8(Cpu_features const & cpu_features,
                        unsigned char * dprofile,
                        unsigned char const * score_matrix,
                        unsigned char const * dseq) -> void;


// One block of cells for the 16-bit width (SSE4.1 unsigned-min path when
// available, else the generic kernel). The 8-bit width has no SSE4.1
// variant, so search8 calls its generic kernels directly.
auto dispatch_align_regular_16(Cpu_features const & cpu_features,
                               VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                               VECTORTYPE16 ** qp,
                               VECTORTYPE16 const * Qm, VECTORTYPE16 const * Rm,
                               uint64_t ql, VECTORTYPE16 const * F0,
                               uint64_t * dir_long, VECTORTYPE16 const * H0) -> void;
auto dispatch_align_masked_16(Cpu_features const & cpu_features,
                              VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                              VECTORTYPE16 ** qp,
                              VECTORTYPE16 const * Qm, VECTORTYPE16 const * Rm,
                              uint64_t ql, VECTORTYPE16 const * F0,
                              uint64_t * dir_long, VECTORTYPE16 const * H0,
                              VECTORTYPE16 const * Mm, VECTORTYPE16 * MQ,
                              VECTORTYPE16 const * MR, VECTORTYPE16 const * MQ0) -> void;


// Generic (SSE2) kernels, defined in search16.cpp / search8.cpp and called
// by the dispatchers above as the always-available fallback.
auto dprofile_fill16(unsigned short * dprofile,
                     unsigned short const * score_matrix,
                     unsigned char const * dseq) -> void;
auto dprofile_fill8(unsigned char * dprofile,
                    unsigned char const * score_matrix,
                    unsigned char const * dseq) -> void;
auto align_cells_regular_16(VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                            VECTORTYPE16 ** qp,
                            VECTORTYPE16 const * Qm, VECTORTYPE16 const * Rm,
                            uint64_t ql, VECTORTYPE16 const * F0,
                            uint64_t * dir_long, VECTORTYPE16 const * H0) -> void;
auto align_cells_masked_16(VECTORTYPE16 * Sm, VECTORTYPE16 * hep,
                           VECTORTYPE16 ** qp,
                           VECTORTYPE16 const * Qm, VECTORTYPE16 const * Rm,
                           uint64_t ql, VECTORTYPE16 const * F0,
                           uint64_t * dir_long, VECTORTYPE16 const * H0,
                           VECTORTYPE16 const * Mm, VECTORTYPE16 * MQ,
                           VECTORTYPE16 const * MR, VECTORTYPE16 const * MQ0) -> void;

#endif
