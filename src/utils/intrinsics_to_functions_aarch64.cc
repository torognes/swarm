/*
  SWARM

  Copyright (C) 2012-2025 Torbjorn Rognes and Frederic Mahe

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

// cancel that, rely on function overloading!! Wait, it is not
// posssible because some x86 functions have the exact same signature
// but use different intrinsics. A solution could be to split x86 into
// two files, one per search8/16. The other arch could use only one
// file, with function overload?

#ifdef __aarch64__

#include <cstdint>  // uint64_t, uint16_t
#include <arm_neon.h>


constexpr uint16x8_t neon_mask16 =
  {0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000};

constexpr uint8x16_t neon_mask8 =
  { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

constexpr uint16x8_t neon_shift16 = { 0, 0, 0, 0, 8, 8, 8, 8 };


auto cast_vector16(uint16_t * ptr) -> uint16_t* {
  // dummy function, needed to match x86-64 code
  return ptr;
}

auto cast_vector8(uint8_t * ptr) -> uint8_t* {
  // dummy function, needed to match x86-64 code
  return ptr;
}

auto cast_vector8(uint8x16_t* ptr) -> uint8x16_t* {
  // dummy function, needed to match x86-64 code
  return ptr;
}

// only used in v_merge_lo8()
auto cast_vector8_real(uint8_t * ptr) -> uint8x16_t* {
  return reinterpret_cast<uint8x16_t*>(ptr);
}

auto cast_vector64(uint8_t * ptr) -> uint64_t* {
  return reinterpret_cast<uint64_t*>(ptr);
}

auto v_load16(uint16_t const * ptr) -> uint16x8_t {
  return vld1q_u16(ptr);
}

// only in search8
auto v_load_64(uint8_t * ptr) -> uint64x2_t {
  return vld1q_dup_u64(cast_vector64(ptr));
}

auto v_store16(uint16_t * ptr, uint16x8_t cpu_register) -> void {
  vst1q_u16(ptr, cpu_register);
}

auto v_store8(uint8_t * ptr, uint8x16_t cpu_register) -> void {
  vst1q_u8(ptr, cpu_register);
}

// only in search8
auto v_merge_lo_8(uint8x16_t lhs, uint8_t& rhs) -> uint8x16_t {
  // vzip1q_u8: interleaves the lower halves of two uint8x16_t
  auto * rhs_ptr = &rhs;
  return vzip1q_u8(lhs, *cast_vector8_real(rhs_ptr));
}

// only in search8
auto v_merge_lo_8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t {
  return vzip1q_u8(lhs, rhs);
}

auto v_merge_lo_16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vzip1q_u16(lhs, rhs);
}

auto v_merge_hi_16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vzip2q_u16(lhs, rhs);
}

auto v_merge_lo_32(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vreinterpretq_u16_u32(
             vzip1q_u32(
                 vreinterpretq_u32_u16(lhs),
                 vreinterpretq_u32_u16(rhs)
             )
         );
}

auto v_merge_hi_32(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vreinterpretq_u16_u32(
             vzip2q_u32(
                 vreinterpretq_u32_u16(lhs),
                 vreinterpretq_u32_u16(rhs)
             )
         );
}

auto v_merge_lo_64(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vreinterpretq_u16_u64(
            vcombine_u64(
                vget_low_u64(vreinterpretq_u64_u16(lhs)),
                vget_low_u64(vreinterpretq_u64_u16(rhs))
            )
         );
}

auto v_merge_hi_64(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vreinterpretq_u16_u64(
             vcombine_u64(
                vget_high_u64(vreinterpretq_u64_u16(lhs)),
                vget_high_u64(vreinterpretq_u64_u16(rhs))
             )
         );
}

auto v_min16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vminq_u16(lhs, rhs);
}

auto v_min8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t {
  return vminq_u8(lhs, rhs);
}

auto v_add16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vqaddq_u16(lhs, rhs);
}

auto v_add8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t {
  return vqaddq_u8(lhs, rhs);
}

auto v_sub16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vqsubq_u16(lhs, rhs);
}

auto v_sub8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t {
  return vqsubq_u8(lhs, rhs);
}
  
auto v_dup16(int16_t value) -> uint16x8_t {
  // broadcast a uint16_t to all elements of destination
  return vdupq_n_u16(static_cast<uint16_t>(value));
}

auto v_dup8(uint8_t value) -> uint8x16_t {
  // broadcast a uint8_t to all elements of destination
  return vdupq_n_u8(value);
}

auto v_zero16() -> uint16x8_t {
  static constexpr int16_t zero = 0;
  return v_dup16(zero);
}

auto v_zero8() -> uint8x16_t {
  static constexpr uint8_t zero = 0;
  return v_dup8(zero);
}

auto v_and16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return vandq_u16(lhs, rhs);
}

auto v_and8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t {
  return vandq_u8(lhs, rhs);
}

auto v_xor16(uint16x8_t lhs, uint16x8_t rhs) -> uint16x8_t {
  return veorq_u16(lhs, rhs);
}

auto v_xor8(uint8x16_t lhs, uint8x16_t rhs) -> uint8x16_t {
  return veorq_u8(lhs, rhs);
}

auto v_shift_left16(uint16x8_t vector) -> uint16x8_t {
  // shift vector to the left by n bytes, pad with zeros
  static constexpr auto n_bytes = 7;
  return vextq_u16(v_zero16(), vector, n_bytes);
}

auto v_shift_left8(uint8x16_t vector) -> uint8x16_t {
  // shift vector to the left by n bytes, pad with zeros
  static constexpr auto n_bytes = 15;
  return vextq_u8(v_zero8(), vector, n_bytes);
}

auto v_mask_eq16(uint16x8_t lhs, uint16x8_t rhs) -> uint16_t {
  // - compare vectors of integers for equality
  // - mask
  // - add and return an uint16
  return vaddvq_u16(vandq_u16(vceqq_u16(lhs, rhs), neon_mask16));
}

auto v_mask_eq8(uint8x16_t lhs, uint8x16_t rhs) -> uint16_t {
  // - compare vectors of integers for equality
  // - Bitwise AND, pairwise add,
  // - shift left by neon_shift16,
  // - add and return an uint16
  return vaddvq_u16(
             vshlq_u16(
                 vpaddlq_u8(
                     vandq_u8(
                         vceqq_u8(lhs, rhs),
                         neon_mask8)),
                 neon_shift16)
         );
}

#endif
