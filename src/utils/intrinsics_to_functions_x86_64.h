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

#ifdef __x86_64__
#ifdef __SSE2__

#include <cstdint>  // uint64_t
#include <emmintrin.h>  // SSE2


auto cast_vector16(uint16_t* ptr) -> __m128i*;

auto cast_vector8(uint8_t* ptr) -> __m128i*;

auto v_load16(uint16_t* ptr) -> __m128i;

auto v_load8(uint8_t* ptr) -> __m128i;

auto v_load_64(uint8_t* ptr) -> __m128i;

auto v_store16(uint16_t* ptr, __m128i cpu_register) -> void;

auto v_store8(uint8_t* ptr, __m128i cpu_register) -> void;

auto v_merge_lo_8(__m128i lhs, __m128i rhs) -> __m128i;

auto v_merge_lo_16(__m128i lhs, __m128i rhs) -> __m128i;

auto v_merge_hi_8(__m128i lhs, __m128i rhs) -> __m128i;

auto v_merge_hi_16(__m128i lhs, __m128i rhs) -> __m128i;

auto v_merge_lo_32(__m128i lhs, __m128i rhs) -> __m128i;

auto v_merge_hi_32(__m128i lhs, __m128i rhs) -> __m128i;

auto v_merge_lo_64(__m128i lhs, __m128i rhs) -> __m128i;

auto v_merge_hi_64(__m128i lhs, __m128i rhs) -> __m128i;

auto v_min16(__m128i lhs, __m128i rhs) -> __m128i;

auto v_min8(__m128i lhs, __m128i rhs) -> __m128i;

auto v_add16(__m128i lhs, __m128i rhs) -> __m128i;

auto v_add8(__m128i lhs, __m128i rhs) -> __m128i;

auto v_sub16(__m128i lhs, __m128i rhs) -> __m128i;

auto v_sub8(__m128i lhs, __m128i rhs) -> __m128i;

auto v_dup16(int16_t value) -> __m128i;

auto v_dup8(int8_t value) -> __m128i;

auto v_zero16() -> __m128i;

auto v_zero8() -> __m128i;

auto v_and16(__m128i lhs, __m128i rhs) -> __m128i;

auto v_and8(__m128i lhs, __m128i rhs) -> __m128i;

auto v_xor16(__m128i lhs, __m128i rhs) -> __m128i;

auto v_xor8(__m128i lhs, __m128i rhs) -> __m128i;

auto v_shift_left16(__m128i vector) -> __m128i;

auto v_shift_left8(__m128i vector) -> __m128i;

auto v_mask_eq16(__m128i lhs, __m128i rhs) -> uint16_t;

auto v_mask_eq8(__m128i lhs, __m128i rhs) -> uint16_t;

#endif
#endif
