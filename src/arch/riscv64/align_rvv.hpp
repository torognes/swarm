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

#ifndef SWARM_ARCH_RISCV64_ALIGN_RVV_H
#define SWARM_ARCH_RISCV64_ALIGN_RVV_H


// Native RVV alignment kernels (align_rvv.cpp), used instead of the
// generic-vector kernel of utils/align_cells.hpp when the V extension is
// enabled. Same arguments, same results, with one difference: the 16-bit
// kernel stores one direction bit per channel (Mask_shape<bits8>), not two.

#ifdef __riscv_vector

#include "intrinsics_to_functions.hpp"  // v_u16_t, v_u8_t
#include <cstdint>  // uint64_t


auto align_rvv_regular_8(v_u8_t * Sm, v_u8_t * hep, v_u8_t ** qp,
                         v_u8_t const & Qm, v_u8_t const & Rm,
                         uint64_t ql, v_u8_t const & F0,
                         uint64_t * dir_long, v_u8_t const & H0) -> void;
auto align_rvv_masked_8(v_u8_t * Sm, v_u8_t * hep, v_u8_t ** qp,
                        v_u8_t const & Qm, v_u8_t const & Rm,
                        uint64_t ql, v_u8_t const & F0,
                        uint64_t * dir_long, v_u8_t const & H0,
                        v_u8_t const * Mm, v_u8_t * MQ,
                        v_u8_t const * MR, v_u8_t const * MQ0) -> void;

auto align_rvv_regular_16(v_u16_t * Sm, v_u16_t * hep, v_u16_t ** qp,
                          v_u16_t const & Qm, v_u16_t const & Rm,
                          uint64_t ql, v_u16_t const & F0,
                          uint64_t * dir_long, v_u16_t const & H0) -> void;
auto align_rvv_masked_16(v_u16_t * Sm, v_u16_t * hep, v_u16_t ** qp,
                         v_u16_t const & Qm, v_u16_t const & Rm,
                         uint64_t ql, v_u16_t const & F0,
                         uint64_t * dir_long, v_u16_t const & H0,
                         v_u16_t const * Mm, v_u16_t * MQ,
                         v_u16_t const * MR, v_u16_t const * MQ0) -> void;

#endif  // __riscv_vector

#endif  // SWARM_ARCH_RISCV64_ALIGN_RVV_H
