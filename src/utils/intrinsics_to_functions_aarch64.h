/*
  SWARM

  Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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

// name functions after input/output type:
// uint64x2_t -> 64
// uint16x8_t -> 16
// uint8x16_t -> 8

#ifdef __aarch64__

#include <cstdint>  // uint64_t, uint16_t
#include <arm_neon.h>


auto cast_vector16(uint16_t * ptr) -> uint16_t*;

auto cast_vector8(uint8_t * ptr) -> uint8_t*;

auto cast_vector8(uint8x16_t* ptr) -> uint8x16_t*;

auto v_load16(uint16_t const * ptr) -> uint16x8_t;

// only in search8
auto v_load_64(uint8_t * ptr) -> uint64x2_t;

auto v_store16(uint16_t * ptr, uint16x8_t cpu_register) -> void;

auto v_store8(uint8_t * ptr, uint8x16_t cpu_register) -> void;

// only in search8
auto v_merge_lo_8(uint8x16_t lhs, uint8_t rhs) -> uint8x16_t;

// only in search8
auto v_merge_lo_8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t;

auto v_merge_lo_16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_merge_hi_16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_merge_lo_32(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_merge_hi_32(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_merge_lo_64(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_merge_hi_64(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_min16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_min8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t;

auto v_add16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_add8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t;

auto v_sub16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_sub8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t;
  
auto v_dup16(uint16_t value) -> uint16x8_t;

auto v_dup8(uint8_t value) -> uint8x16_t;

auto v_zero16() -> uint16x8_t;

auto v_zero8() -> uint8x16_t;

auto v_and16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_and8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t;

auto v_xor16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t;

auto v_xor8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t;

auto v_shift_left16(uint16x8_t vector) -> uint16x8_t;

auto v_shift_left8(uint8x16_t vector) -> uint8x16_t;

auto v_mask_eq16(uint16x8_t lhs, uint16x8_t rhs) -> uint16_t;

auto v_mask_eq8(uint8x16_t lhs, uint8x16_t rhs) -> uint16_t;

#endif
