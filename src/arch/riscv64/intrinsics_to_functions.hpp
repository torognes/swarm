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

#ifndef SWARM_ARCH_RISCV64_INTRINSICS_TO_FUNCTIONS_H
#define SWARM_ARCH_RISCV64_INTRINSICS_TO_FUNCTIONS_H


// name functions after input/output type:
// v_u16_t -> 16
// v_u8_t -> 8

// RISC-V: no intrinsics here. The RVV intrinsic types are sizeless
// (their length is only known at run time), so they cannot be stored
// in the VECTORTYPE arrays the search kernels use. These functions are
// written instead with GCC/clang generic vectors (vector_size), fixed
// at 128 bits like the other ports; with -march=rv64gcv the compiler
// lowers them to RVV instructions (any VLEN >= 128). Nothing below is
// specific to RISC-V.

#include <cstdint>  // uint16_t, uint8_t, int16_t


using v_u16_t __attribute__((vector_size(16))) = uint16_t;
using v_u8_t __attribute__((vector_size(16))) = uint8_t;


auto cast_vector16(uint16_t * ptr) -> uint16_t*;
auto cast_vector16(uint16_t const * ptr) -> uint16_t const *;

auto cast_vector8(uint8_t * ptr) -> uint8_t*;
auto cast_vector8(uint8_t const * ptr) -> uint8_t const *;

auto v_load16(uint16_t const * ptr) -> v_u16_t;

// only in search8
auto v_load_64(uint8_t const * ptr) -> v_u8_t;

auto v_store16(uint16_t * ptr, v_u16_t cpu_register) -> void;

auto v_store8(uint8_t * ptr, v_u8_t cpu_register) -> void;

// only in search8
auto v_merge_lo_8(v_u8_t lhs, uint8_t const & rhs) -> v_u8_t;

auto v_merge_lo_16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_merge_lo_16(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_merge_hi_16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_merge_hi_16(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_merge_lo_32(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_merge_lo_32(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_merge_hi_32(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_merge_hi_32(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_merge_lo_64(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_merge_lo_64(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_merge_hi_64(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_merge_hi_64(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_min16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_min8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_add16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_add8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_sub16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_sub8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_dup16(int16_t value) -> v_u16_t;

auto v_dup8(uint8_t value) -> v_u8_t;

auto v_zero16() -> v_u16_t;

auto v_zero8() -> v_u8_t;

auto v_and16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_and8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_xor16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t;

auto v_xor8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t;

auto v_shift_left16(v_u16_t vector) -> v_u16_t;

auto v_shift_left8(v_u8_t vector) -> v_u8_t;

auto v_mask_eq16(v_u16_t lhs, v_u16_t rhs) -> uint16_t;

auto v_mask_eq8(v_u8_t lhs, v_u8_t rhs) -> uint16_t;

#endif  // SWARM_ARCH_RISCV64_INTRINSICS_TO_FUNCTIONS_H
