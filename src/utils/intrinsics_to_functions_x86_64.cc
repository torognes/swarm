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

#ifdef __x86_64__
#ifdef __SSE2__

#include <cstdint>  // uint64_t
#include <emmintrin.h>  // SSE2


auto cast_vector16(uint16_t* ptr) -> __m128i* {
  return reinterpret_cast<__m128i*>(ptr);
}

auto cast_vector8(uint8_t* ptr) -> __m128i* {
  return reinterpret_cast<__m128i*>(ptr);
}

auto v_load16(uint16_t* ptr) -> __m128i {
  return _mm_load_si128(cast_vector16(ptr));
}

// only in ssse3
auto v_load8(uint8_t* ptr) -> __m128i {
  return _mm_load_si128(cast_vector8(ptr));
}

// only in search8
auto v_load_64(uint8_t* ptr) -> __m128i {
  return _mm_loadl_epi64(cast_vector8(ptr));
}

auto v_store16(uint16_t* ptr, __m128i cpu_register) -> void {
  _mm_store_si128(cast_vector16(ptr), cpu_register);
}

auto v_store8(uint8_t* ptr, __m128i cpu_register) -> void {
  _mm_store_si128(cast_vector8(ptr), cpu_register);
}

// only in search8
auto v_merge_lo_8(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_unpacklo_epi8(lhs, rhs);
}

auto v_merge_lo_16(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_unpacklo_epi16(lhs, rhs);
}

auto v_merge_hi_16(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_unpackhi_epi16(lhs, rhs);
}

auto v_merge_lo_32(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_unpacklo_epi32(lhs, rhs);
}

auto v_merge_hi_32(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_unpackhi_epi32(lhs, rhs);
}

auto v_merge_lo_64(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_unpacklo_epi64(lhs, rhs);
}

auto v_merge_hi_64(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_unpackhi_epi64(lhs, rhs);
}

// refactoring: why more complicated than vmin8()? _mm_min_epu16??
auto v_min16(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_subs_epu16(lhs, _mm_subs_epu16(lhs, rhs));
}

auto v_min8(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_min_epu8(lhs, rhs);
}

auto v_add16(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_adds_epu16(lhs, rhs);
}

auto v_add8(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_adds_epu8(lhs, rhs);
}

auto v_sub16(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_subs_epu16(lhs, rhs);
}

auto v_sub8(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_subs_epu8(lhs, rhs);
}

auto v_dup16(int16_t value) -> __m128i {
  // broadcast a 16-bit integer to all elements of destination
  return _mm_set1_epi16(value); // takes a short int (int16_t)
}

auto v_dup8(int8_t value) -> __m128i {
  // broadcast a 8-bit integer to all elements of destination
  return _mm_set1_epi8(value); // takes a char (int8_t)
}

auto v_zero16() -> __m128i {
  static constexpr int16_t zero = 0;
  return v_dup16(zero);
}

auto v_zero8() -> __m128i {
  static constexpr int8_t zero = 0;
  return v_dup8(zero);
}

auto v_and16(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_and_si128(lhs, rhs);
}

auto v_and8(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_and_si128(lhs, rhs);
}

auto v_xor16(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_xor_si128(lhs, rhs);
}

auto v_xor8(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_xor_si128(lhs, rhs);
}

auto v_shift_left16(__m128i vector) -> __m128i {
  // shift vector to the left by n bytes, pad with zeros
  static constexpr auto n_bytes = 2;
  return _mm_slli_si128(vector, n_bytes);
}

auto v_shift_left8(__m128i vector) -> __m128i {
  // shift vector to the left by n bytes, pad with zeros
  static constexpr auto n_bytes = 1;
  return _mm_slli_si128(vector, n_bytes);
}

auto v_mask_eq16(__m128i lhs, __m128i rhs) -> uint16_t {
  // - compare vectors of 16-bit integers for equality
  // - create mask from the most significant bit of each 8-bit
  //   element,
  // - return an int
  return static_cast<uint16_t>(_mm_movemask_epi8(_mm_cmpeq_epi16(lhs, rhs)));
}

auto v_mask_eq8(__m128i lhs, __m128i rhs) -> uint16_t {
  // - compare vectors of 8-bit integers for equality
  // - create mask from the most significant bit of each 8-bit
  //   element,
  // - return an int
  return static_cast<uint16_t>(_mm_movemask_epi8(_mm_cmpeq_epi8(lhs, rhs)));
}

#endif
#endif
