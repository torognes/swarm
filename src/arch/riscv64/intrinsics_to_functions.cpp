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

// name functions after input/output type:
// v_u16_t -> 16
// v_u8_t -> 8

// Generic-vector implementation (see intrinsics_to_functions.hpp).
// Lanes are addressed by index, as in the x86_64 SSE2 reference: lane
// 0 is the lowest address. Shuffle indices 0-15 (or 0-7) select from
// lhs, 16-31 (or 8-15) from rhs. Byte-vector merges at a wider
// granularity shuffle byte pairs, quads or octets, so they do not
// depend on the byte order.

#if !defined(__clang__) && defined(__GNUC__) && (__GNUC__ < 12)
#error GCC 12 or more recent is required (__builtin_shufflevector)
#endif

#include "intrinsics_to_functions.hpp"
#include "../../utils/simd_alignment.hpp"  // simd_vector_bytes
#include <array>
#include <cstdint>  // uint16_t, uint64_t, uint8_t, int16_t
#include <cstring>  // std::memcpy


static_assert(sizeof(v_u8_t) == simd_vector_bytes, "v_u8_t must be 128-bit");
static_assert(sizeof(v_u16_t) == simd_vector_bytes, "v_u16_t must be 128-bit");


namespace {

  // lane counts from simd_vector_bytes rather than sizeof(v_u8_t):
  // cppcheck does not know vector_size types and evaluates their
  // sizeof as the size of one element
  constexpr auto n_lanes16 = simd_vector_bytes / sizeof(uint16_t);  // 8
  constexpr auto n_half_lanes8 = simd_vector_bytes / 2;  // 8

  // lanes of a comparison result are all-ones or all-zeros; reuse
  // them as an unsigned mask. Template: the result of comparing two
  // generic vectors is a signed vector whose exact element type
  // differs between compilers.
  template <typename Comparison>
  auto mask8(Comparison const comparison) -> v_u8_t {
    static_assert(sizeof(Comparison) == sizeof(v_u8_t), "not a 128-bit vector");
    return reinterpret_cast<v_u8_t>(comparison);
  }

  template <typename Comparison>
  auto mask16(Comparison const comparison) -> v_u16_t {
    static_assert(sizeof(Comparison) == sizeof(v_u16_t), "not a 128-bit vector");
    return reinterpret_cast<v_u16_t>(comparison);
  }

  // movemask: weight of each lane in the bitmask returned by
  // v_mask_eq8 (one bit per byte) and v_mask_eq16 (two bits per
  // 16-bit lane), as in _mm_movemask_epi8()
  constexpr v_u8_t lane_bits8 = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
                                  0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

  constexpr v_u16_t lane_bits16 = { 0x0003, 0x000c, 0x0030, 0x00c0,
                                    0x0300, 0x0c00, 0x3000, 0xc000 };

}  // end of anonymous namespace


auto cast_vector16(uint16_t * const ptr) -> uint16_t* {
  // dummy function, needed to match x86-64 code
  return ptr;
}

auto cast_vector16(uint16_t const * const ptr) -> uint16_t const * {
  // dummy function, needed to match x86-64 code
  return ptr;
}

auto cast_vector8(uint8_t * const ptr) -> uint8_t* {
  // dummy function, needed to match x86-64 code
  return ptr;
}

auto cast_vector8(uint8_t const * const ptr) -> uint8_t const * {
  // dummy function, needed to match x86-64 code
  return ptr;
}

// std::memcpy: loads and stores go through bytes, so they need neither
// alignment nor type punning; the compiler folds them into single
// vector loads and stores

auto v_load16(uint16_t const * const ptr) -> v_u16_t {
  v_u16_t result;
  std::memcpy(&result, ptr, sizeof(result));
  return result;
}

// only in search8
auto v_load_64(uint8_t const * const ptr) -> v_u8_t {
  // load 64 bits into both halves of the vector
  uint64_t half {0};
  std::memcpy(&half, ptr, sizeof(half));
  std::array<uint64_t, 2> const halves {{half, half}};
  v_u8_t result;
  std::memcpy(&result, halves.data(), sizeof(result));
  return result;
}

auto v_store16(uint16_t * const ptr, v_u16_t const cpu_register) -> void {
  std::memcpy(ptr, &cpu_register, sizeof(cpu_register));
}

auto v_store8(uint8_t * const ptr, v_u8_t const cpu_register) -> void {
  std::memcpy(ptr, &cpu_register, sizeof(cpu_register));
}

// only in search8
auto v_merge_lo_8(v_u8_t const lhs, uint8_t const & rhs) -> v_u8_t {
  // interleave the low halves of lhs and of the 64 bits at &rhs (rhs
  // is the first of 8 bytes, not a lone byte; the shuffle ignores the
  // duplicated high half)
  return __builtin_shufflevector(lhs, v_load_64(&rhs),
                                 0, 16, 1, 17, 2, 18, 3, 19,
                                 4, 20, 5, 21, 6, 22, 7, 23);
}

auto v_merge_lo_16(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return __builtin_shufflevector(lhs, rhs, 0, 8, 1, 9, 2, 10, 3, 11);
}

auto v_merge_lo_16(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return __builtin_shufflevector(lhs, rhs,
                                 0, 1, 16, 17, 2, 3, 18, 19,
                                 4, 5, 20, 21, 6, 7, 22, 23);
}

auto v_merge_hi_16(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return __builtin_shufflevector(lhs, rhs, 4, 12, 5, 13, 6, 14, 7, 15);
}

auto v_merge_hi_16(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return __builtin_shufflevector(lhs, rhs,
                                 8, 9, 24, 25, 10, 11, 26, 27,
                                 12, 13, 28, 29, 14, 15, 30, 31);
}

auto v_merge_lo_32(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return __builtin_shufflevector(lhs, rhs, 0, 1, 8, 9, 2, 3, 10, 11);
}

auto v_merge_lo_32(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return __builtin_shufflevector(lhs, rhs,
                                 0, 1, 2, 3, 16, 17, 18, 19,
                                 4, 5, 6, 7, 20, 21, 22, 23);
}

auto v_merge_hi_32(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return __builtin_shufflevector(lhs, rhs, 4, 5, 12, 13, 6, 7, 14, 15);
}

auto v_merge_hi_32(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return __builtin_shufflevector(lhs, rhs,
                                 8, 9, 10, 11, 24, 25, 26, 27,
                                 12, 13, 14, 15, 28, 29, 30, 31);
}

auto v_merge_lo_64(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return __builtin_shufflevector(lhs, rhs, 0, 1, 2, 3, 8, 9, 10, 11);
}

auto v_merge_lo_64(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return __builtin_shufflevector(lhs, rhs,
                                 0, 1, 2, 3, 4, 5, 6, 7,
                                 16, 17, 18, 19, 20, 21, 22, 23);
}

auto v_merge_hi_64(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return __builtin_shufflevector(lhs, rhs, 4, 5, 6, 7, 12, 13, 14, 15);
}

auto v_merge_hi_64(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return __builtin_shufflevector(lhs, rhs,
                                 8, 9, 10, 11, 12, 13, 14, 15,
                                 24, 25, 26, 27, 28, 29, 30, 31);
}

// min, saturating add and saturating sub (unsigned): written with
// comparison masks, which GCC and clang both accept on generic
// vectors in C++11 (the vector ?: operator is not portable)

auto v_min16(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return rhs ^ ((lhs ^ rhs) & mask16(lhs < rhs));
}

auto v_min8(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return rhs ^ ((lhs ^ rhs) & mask8(lhs < rhs));
}

auto v_add16(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  // a lane that wraps around becomes smaller than lhs: clamp it to max
  auto const sum = lhs + rhs;
  return sum | mask16(sum < lhs);
}

auto v_add8(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  // a lane that wraps around becomes smaller than lhs: clamp it to max
  auto const sum = lhs + rhs;
  return sum | mask8(sum < lhs);
}

auto v_sub16(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  // a lane that would go below zero is clamped to zero
  return (lhs - rhs) & mask16(lhs > rhs);
}

auto v_sub8(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  // a lane that would go below zero is clamped to zero
  return (lhs - rhs) & mask8(lhs > rhs);
}

auto v_dup16(int16_t const value) -> v_u16_t {
  // broadcast a uint16_t to all elements of destination
  auto const lane = static_cast<uint16_t>(value);
  v_u16_t const result = { lane, lane, lane, lane, lane, lane, lane, lane };
  return result;
}

auto v_dup8(uint8_t const value) -> v_u8_t {
  // broadcast a uint8_t to all elements of destination
  v_u8_t const result = { value, value, value, value, value, value, value, value,
                          value, value, value, value, value, value, value, value };
  return result;
}

auto v_zero16() -> v_u16_t {
  static constexpr int16_t zero = 0;
  return v_dup16(zero);
}

auto v_zero8() -> v_u8_t {
  static constexpr uint8_t zero = 0;
  return v_dup8(zero);
}

auto v_and16(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return lhs & rhs;
}

auto v_and8(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return lhs & rhs;
}

auto v_xor16(v_u16_t const lhs, v_u16_t const rhs) -> v_u16_t {
  return lhs ^ rhs;
}

auto v_xor8(v_u8_t const lhs, v_u8_t const rhs) -> v_u8_t {
  return lhs ^ rhs;
}

auto v_shift_left16(v_u16_t const vector) -> v_u16_t {
  // shift vector to the left by one lane (towards higher addresses),
  // pad with zeros
  return __builtin_shufflevector(vector, v_zero16(), 8, 0, 1, 2, 3, 4, 5, 6);
}

auto v_shift_left8(v_u8_t const vector) -> v_u8_t {
  // shift vector to the left by one lane (towards higher addresses),
  // pad with zeros
  return __builtin_shufflevector(vector, v_zero8(),
                                 16, 0, 1, 2, 3, 4, 5, 6,
                                 7, 8, 9, 10, 11, 12, 13, 14);
}

auto v_mask_eq16(v_u16_t const lhs, v_u16_t const rhs) -> uint16_t {
  // - compare vectors of integers for equality
  // - keep the two bits of each lane's weight
  // - combine the lanes and return an uint16
  auto const bits = mask16(lhs == rhs) & lane_bits16;
  uint16_t result {0};
  for (auto i = 0U; i < n_lanes16; ++i) {
    result = static_cast<uint16_t>(result | bits[i]);
  }
  return result;
}

auto v_mask_eq8(v_u8_t const lhs, v_u8_t const rhs) -> uint16_t {
  // - compare vectors of integers for equality
  // - keep each lane's bit
  // - combine the low and high halves and return an uint16
  static constexpr auto bits_per_byte = 8U;
  auto const bits = mask8(lhs == rhs) & lane_bits8;
  unsigned int low_half {0};
  unsigned int high_half {0};
  for (auto i = 0U; i < n_half_lanes8; ++i) {
    low_half |= bits[i];
    high_half |= bits[i + n_half_lanes8];
  }
  return static_cast<uint16_t>(low_half | (high_half << bits_per_byte));
}
