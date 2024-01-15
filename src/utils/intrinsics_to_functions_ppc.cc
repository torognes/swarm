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

#ifdef __PPC__

#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#include <cstdint>  // uint64_t, uint16_t

using v_u64_t = vector unsigned long long;
using v_u32_t = vector unsigned int;
using v_u16_t = vector unsigned short;
using v_u8_t = vector unsigned char;

constexpr v_u8_t perm_merge_long_low =
  {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
   0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17};

constexpr v_u8_t perm_merge_long_high =
  {0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
   0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

constexpr v_u8_t perm_bits =
  {0x78, 0x70, 0x68, 0x60, 0x58, 0x50, 0x48, 0x40,
   0x38, 0x30, 0x28, 0x20, 0x18, 0x10, 0x08, 0x00};

// source: power vector intrinsic programming reference (chapter 4)

auto cast_vector16(uint16_t * ptr) -> v_u16_t* {
  return reinterpret_cast<v_u16_t*>(ptr);
}

auto cast_vector8(uint8_t * ptr) -> v_u8_t* {
  return reinterpret_cast<v_u8_t*>(ptr);
}

auto v_load16(uint16_t const * ptr) -> v_u16_t {
  return vec_splats(*ptr);  // dereference
}

// only in search8
auto v_load_64(uint64_t const * ptr) -> v_u64_t {
  return vec_splats(*ptr);
}

auto v_store16(uint16_t * ptr, v_u16_t cpu_register) -> void {
  // store a vector at ptr address + displacement
  static constexpr auto displacement = 0LL;
  // useless in the original macro?
  // - (VECTORTYPE)(cpu_register)
  // - CAST_VECTOR_p(c)
  vec_st(cpu_register, displacement, ptr);
}

auto v_store8(uint8_t * ptr, v_u8_t cpu_register) -> void {
  // store a vector at ptr address + displacement
  static constexpr auto displacement = 0LL;
  // useless in the original macro?
  // - (VECTORTYPE)(cpu_register)
  // - CAST_VECTOR_p(c)
  vec_st(cpu_register, displacement, ptr);
}

// only in search8
auto v_merge_lo_8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  return vec_mergeh(lhs, rhs);
}

auto v_merge_lo_16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  return vec_mergeh(lhs, rhs);
}

auto v_merge_lo_16(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  return vec_mergeh(lhs, rhs);
}

auto v_merge_hi_16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  return vec_mergel(lhs, rhs);
}

auto v_merge_hi_16(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  // search8: v_merge_hi_16(a, b) (VECTORTYPE)vec_mergel((vector short)(a), (vector short)(b))
  // refactoring: try without casting
  return static_cast<v_u8_t>(
             vec_mergel(static_cast<v_u16_t>(lhs),
                        static_cast<v_u16_t>(rhs))
         );
}

auto v_merge_lo_32(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  // search16: v_merge_lo_32(a, b) (VECTORTYPE)vec_mergeh((vector int)(a), (vector int)(b))
  // decision: was casting to signed int; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u16_t>(
             vec_mergeh(static_cast<v_u32_t>(lhs),
                        static_cast<v_u32_t>(rhs))
         );
}

auto v_merge_lo_32(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  // search8: v_merge_lo_32(a, b) (VECTORTYPE)vec_mergeh((vector int)(a), (vector int)(b))
  // decision: was casting to signed int; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u8_t>(
             vec_mergeh(static_cast<v_u32_t>(lhs),
                        static_cast<v_u32_t>(rhs))
         );
}

auto v_merge_hi_32(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  // search16: v_merge_hi_32(a, b) (VECTORTYPE)vec_mergel((vector int)(a), (vector int)(b))
  // decision: was casting to signed int; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u16_t>(
             vec_mergel(static_cast<v_u32_t>(lhs),
                        static_cast<v_u32_t>(rhs))
         );
}

auto v_merge_hi_32(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  // search8: v_merge_hi_32(a, b) (VECTORTYPE)vec_mergel((vector int)(a), (vector int)(b))
  // decision: was casting to signed int; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u8_t>(
             vec_mergel(static_cast<v_u32_t>(lhs),
                        static_cast<v_u32_t>(rhs))
         );
}

auto v_merge_lo_64(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  // search16: v_merge_lo_64(a, b) (VECTORTYPE)vec_perm((vector long long)(a), (vector long long)(b), perm_merge_long_low)
  // decision: was casting to signed long long; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u16_t>(
             vec_perm(static_cast<v_u64_t>(lhs),
                      static_cast<v_u64_t>(rhs),
                      perm_merge_long_low)
         );
}

auto v_merge_lo_64(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  // search8: v_merge_lo_64(a, b) (VECTORTYPE)vec_perm((vector long long)(a), (vector long long)(b), perm_merge_long_low)
  // decision: was casting to signed long long; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u8_t>(
             vec_perm(static_cast<v_u64_t>(lhs),
                      static_cast<v_u64_t>(rhs),
                      perm_merge_long_low)
         );
}


auto v_merge_hi_64(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  // search16: v_merge_hi_64(a, b) (VECTORTYPE)vec_perm((vector long long)(a), (vector long long)(b), perm_merge_long_high)
  // decision: was casting to signed long long; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u16_t>(
             vec_perm(static_cast<v_u64_t>(lhs),
                      static_cast<v_u64_t>(rhs),
                      perm_merge_long_high)
         );
}

auto v_merge_hi_64(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  // search8: v_merge_hi_64(a, b) (VECTORTYPE)vec_perm((vector long long)(a), (vector long long)(b), perm_merge_long_high)
  // decision: was casting to signed long long; I use unsigned
  // refactoring: try without casting
  return static_cast<v_u8_t>(
             vec_perm(static_cast<v_u64_t>(lhs),
                      static_cast<v_u64_t>(rhs),
                      perm_merge_long_high)
         );
}

auto v_min16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  return vec_min(lhs, rhs);
}

auto v_min8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  return vec_min(lhs, rhs);
}

auto v_add16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  return vec_adds(lhs, rhs);
}

auto v_add8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  return vec_adds(lhs, rhs);
}

auto v_sub16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  return vec_subs(lhs, rhs);
}

auto v_sub8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  return vec_subs(lhs, rhs);
}

auto v_dup16(uint16_t value) -> v_u16_t {
  // broadcast a uint16_t to all elements of destination
  return vec_splats(value);
}

auto v_dup8(uint8_t value) -> v_u8_t {
  // broadcast a uint8_t to all elements of destination
  return vec_splats(value);
}

auto v_zero16() -> v_u16_t {
  // decision: was using vec_splat_u16, I use v_dup16()
  static constexpr uint16_t zero = 0;
  return v_dup16(zero);
}

auto v_zero8() -> v_u8_t {
  // decision: was using vec_splat_u8, I use v_dup8()
  static constexpr uint8_t zero = 0;
  return v_dup8(zero);
}

auto v_and16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  return vec_and(lhs, rhs);
}

auto v_and8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  return vec_and(lhs, rhs);
}

auto v_xor16(v_u16_t lhs, v_u16_t rhs) -> v_u16_t {
  return vec_xor(lhs, rhs);
}

auto v_xor8(v_u8_t lhs, v_u8_t rhs) -> v_u8_t {
  return vec_xor(lhs, rhs);
}

auto v_shift_left16(v_u16_t vector) -> v_u16_t {
  // shift vector to the left by n bytes, pad with zeros
  // n: 4-bit unsigned literal (in the range 0–15)
  // static constexpr auto n_bytes = 2;
  return vec_sld(vector, v_zero16(), 2);
}

auto v_shift_left8(v_u8_t vector) -> v_u8_t {
  // shift vector to the left by n bytes, pad with zeros
  // n: 4-bit unsigned literal (in the range 0–15)
  // static constexpr auto n_bytes = 1;
  return vec_sld(vector, v_zero8(), 1);
}

auto v_mask_eq16(v_u8_t lhs, v_u8_t rhs) -> uint16_t {
  // - compare vectors of integers for equality
  // - vec_cmpeq -> vector bool char
  // - permute (vec_bperm -> vector unsigned char)
  // - extract a uint16!?
  // refactoring: should function return a uint8_t?
  // refactoring: vec_vbpermq is vec_bperm in the
  // documentation. Modify?
  static constexpr auto fifth_item = 4U;
  return static_cast<v_u16_t>(
             vec_vbpermq(
                 static_cast<v_u8_t>(
                     vec_cmpeq(lhs, rhs)
                     ), perm_bits)
         )[fifth_item];
}
    
auto v_mask_eq8(v_u8_t lhs, v_u8_t rhs) -> uint16_t {
  // - compare vectors of integers for equality
  // - vec_cmpeq -> vector bool char
  // - permute (vec_bperm -> vector unsigned char)
  // - extract a uint16!?
  // refactoring: should function return a uint8_t?
  // refactoring: vec_vbpermq is vec_bperm in the
  // documentation. Modify?
  static constexpr auto fifth_item = 4U;
  return static_cast<v_u16_t>(
             vec_vbpermq(
                 static_cast<v_u8_t>(
                     vec_cmpeq(lhs, rhs)
                     ), perm_bits)
         )[fifth_item];
}

#endif
