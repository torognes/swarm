/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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

#include "swarm.h"
#include "db.h"
#include <array>

// refactoring: C++26 std::simd
#ifdef __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#endif

#include "ssse3.h"

#endif

#include "utils/nt_codec.h"
#include <cstdint>  // int64_t, uint64_t


#ifdef __aarch64__
#include <arm_neon.h>
#elif defined __x86_64__

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#else

#error Unknown architecture
#endif


constexpr unsigned int channels {16};
constexpr unsigned int cdepth {4};
using BYTE = unsigned char;


/* uses 16 unsigned 8-bit values */

#ifdef __aarch64__

using VECTORTYPE = int8x16_t;

#define CAST_VECTOR_p(x) (reinterpret_cast<VECTORTYPE *>(x))

const uint8x16_t neon_mask =
  { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80,
    0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80 };

const uint16x8_t neon_shift = { 0, 0, 0, 0, 8, 8, 8, 8 };

#define v_load_64(a) vld1q_dup_u64((const uint64_t *)(a))
#define v_store(a, b) vst1q_s8((int8_t *)(a), (b))
#define v_merge_lo_8(a, b) vzip1q_s8((a),(b))
#define v_merge_lo_16(a, b) vzip1q_s16((a),(b))
#define v_merge_hi_16(a, b) vzip2q_s16((a),(b))
#define v_merge_lo_32(a, b) vreinterpretq_s16_s32(vzip1q_s32 \
          (vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)))
#define v_merge_hi_32(a, b) vreinterpretq_s16_s32(vzip2q_s32 \
          (vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)))
#define v_merge_lo_64(a, b) vreinterpretq_s16_s64(vcombine_s64 \
          (vget_low_s64(vreinterpretq_s64_s16(a)), \
           vget_low_s64(vreinterpretq_s64_s16(b))))
#define v_merge_hi_64(a, b) vreinterpretq_s16_s64(vcombine_s64 \
          (vget_high_s64(vreinterpretq_s64_s16(a)), \
           vget_high_s64(vreinterpretq_s64_s16(b))))
#define v_min(a, b) vminq_u8((a), (b))
#define v_add(a, b) vqaddq_u8((a), (b))
#define v_sub(a, b) vqsubq_u8((a), (b))
#define v_dup(a) vdupq_n_u8(a)
#define v_zero v_dup(0)
#define v_and(a, b) vandq_u8((a), (b))
#define v_xor(a, b) veorq_u8((a), (b))
#define v_shift_left(a) vextq_u8((v_zero), (a), 15)
#define v_mask_eq(a, b) vaddvq_u16(vshlq_u16(vpaddlq_u8(vandq_u8 \
          ((vceqq_s8((a), (b))), neon_mask)), neon_shift))

#elif defined __x86_64__

using VECTORTYPE = __m128i;

#define CAST_VECTOR_p(x) (reinterpret_cast<VECTORTYPE *>(x))

#define v_load_64(a) _mm_loadl_epi64(CAST_VECTOR_p(a))
#define v_store(a, b) _mm_store_si128(CAST_VECTOR_p(a), (b))
#define v_merge_lo_8(a, b) _mm_unpacklo_epi8((a),(b))
#define v_merge_lo_16(a, b) _mm_unpacklo_epi16((a),(b))
#define v_merge_hi_16(a, b) _mm_unpackhi_epi16((a),(b))
#define v_merge_lo_32(a, b) _mm_unpacklo_epi32((a),(b))
#define v_merge_hi_32(a, b) _mm_unpackhi_epi32((a),(b))
#define v_merge_lo_64(a, b) _mm_unpacklo_epi64((a),(b))
#define v_merge_hi_64(a, b) _mm_unpackhi_epi64((a),(b))
#define v_min(a, b) _mm_min_epu8((a), (b))
#define v_add(a, b) _mm_adds_epu8((a), (b))
#define v_sub(a, b) _mm_subs_epu8((a), (b))
#define v_dup(a) _mm_set1_epi8(a)
#define v_zero v_dup(0)
#define v_and(a, b) _mm_and_si128((a), (b))
#define v_xor(a, b) _mm_xor_si128((a), (b))
#define v_shift_left(a) _mm_slli_si128((a), 1)
#define v_mask_eq(a, b) static_cast<unsigned short>(_mm_movemask_epi8(_mm_cmpeq_epi8((a), (b))))

#elif defined __PPC__

using VECTORTYPE = vector unsigned char;

#define CAST_VECTOR_p(x) reinterpret_cast<VECTORTYPE *>(x)

const vector unsigned char perm_merge_long_low =
  {0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
   0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17};

const vector unsigned char perm_merge_long_high =
  {0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
   0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f};

const vector unsigned char perm_bits =
  { 0x78, 0x70, 0x68, 0x60, 0x58, 0x50, 0x48, 0x40,
    0x38, 0x30, 0x28, 0x20, 0x18, 0x10, 0x08, 0x00 };

#define v_load_64(a) (VECTORTYPE)vec_splats(*((unsigned long long *)(a)))
#define v_store(a, b) vec_st((VECTORTYPE)(b), 0, (VECTORTYPE*)(a))
#define v_merge_lo_8(a, b) vec_mergeh((a), (b))
#define v_merge_lo_16(a, b) (VECTORTYPE)vec_mergeh((vector short)(a),\
                                                   (vector short)(b))
#define v_merge_hi_16(a, b) (VECTORTYPE)vec_mergel((vector short)(a),\
                                                   (vector short)(b))
#define v_merge_lo_32(a, b) (VECTORTYPE)vec_mergeh((vector int)(a), \
                                                   (vector int)(b))
#define v_merge_hi_32(a, b) (VECTORTYPE)vec_mergel((vector int)(a), \
                                                   (vector int)(b))
#define v_merge_lo_64(a, b) (VECTORTYPE)vec_perm((vector long long)(a), \
                                                 (vector long long)(b), \
                                                 perm_merge_long_low)
#define v_merge_hi_64(a, b) (VECTORTYPE)vec_perm((vector long long)(a), \
                                                 (vector long long)(b), \
                                                 perm_merge_long_high)
#define v_min(a, b) vec_min((a), (b))
#define v_add(a, b) vec_adds((a), (b))
#define v_sub(a, b) vec_subs((a), (b))
#define v_dup(a) vec_splats((unsigned char)(a));
#define v_zero vec_splat_u8(0)
#define v_and(a, b) vec_and((a), (b))
#define v_xor(a, b) vec_xor((a), (b))
#define v_shift_left(a) vec_sld((a), v_zero, 1)
#define v_mask_eq(a, b) ((vector unsigned short) \
                         vec_vbpermq((vector unsigned char)             \
                                     vec_cmpeq((a), (b)), perm_bits))[4]

#else

#error Unknown Architecture

#endif


inline void dprofile_fill8(BYTE * dprofile,
                           BYTE * score_matrix,
                           BYTE * dseq)
{
  static constexpr unsigned int multiplier {5};
  VECTORTYPE reg0;
  VECTORTYPE reg1;
  VECTORTYPE reg2;
  VECTORTYPE reg3;
  VECTORTYPE reg4;
  VECTORTYPE reg5;
  VECTORTYPE reg6;
  VECTORTYPE reg7;
  VECTORTYPE reg8;
  VECTORTYPE reg9;
  VECTORTYPE reg10;
  VECTORTYPE reg11;
  VECTORTYPE reg12;
  VECTORTYPE reg13;
  VECTORTYPE reg14;
  VECTORTYPE reg15;

  for(auto j = 0U; j < cdepth; j++)
    {
      std::array<unsigned int, channels> d {{}};
      for(auto i = 0U; i < channels; i++) {
        d[i] = (static_cast<unsigned int>(dseq[j * channels + i])) << multiplier;
      }

      reg0  = v_load_64(score_matrix + d[ 0]);
      reg2  = v_load_64(score_matrix + d[ 2]);
      reg4  = v_load_64(score_matrix + d[ 4]);
      reg6  = v_load_64(score_matrix + d[ 6]);
      reg8  = v_load_64(score_matrix + d[ 8]);
      reg10 = v_load_64(score_matrix + d[10]);
      reg12 = v_load_64(score_matrix + d[12]);
      reg14 = v_load_64(score_matrix + d[14]);

      reg0  = v_merge_lo_8(reg0,  *CAST_VECTOR_p(score_matrix + d[ 1]));
      reg2  = v_merge_lo_8(reg2,  *CAST_VECTOR_p(score_matrix + d[ 3]));
      reg4  = v_merge_lo_8(reg4,  *CAST_VECTOR_p(score_matrix + d[ 5]));
      reg6  = v_merge_lo_8(reg6,  *CAST_VECTOR_p(score_matrix + d[ 7]));
      reg8  = v_merge_lo_8(reg8,  *CAST_VECTOR_p(score_matrix + d[ 9]));
      reg10 = v_merge_lo_8(reg10, *CAST_VECTOR_p(score_matrix + d[11]));
      reg12 = v_merge_lo_8(reg12, *CAST_VECTOR_p(score_matrix + d[13]));
      reg14 = v_merge_lo_8(reg14, *CAST_VECTOR_p(score_matrix + d[15]));

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile + 16 * j +   0, reg0);
      v_store(dprofile + 16 * j +  64, reg3);
      v_store(dprofile + 16 * j + 128, reg2);
      v_store(dprofile + 16 * j + 192, reg7);
      v_store(dprofile + 16 * j + 256, reg1);
      v_store(dprofile + 16 * j + 320, reg11);
      v_store(dprofile + 16 * j + 384, reg6);
      v_store(dprofile + 16 * j + 448, reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(score_matrix + 8 + d[0 ]);
      reg1  = v_load_64(score_matrix + 8 + d[1 ]);
      reg2  = v_load_64(score_matrix + 8 + d[2 ]);
      reg3  = v_load_64(score_matrix + 8 + d[3 ]);
      reg4  = v_load_64(score_matrix + 8 + d[4 ]);
      reg5  = v_load_64(score_matrix + 8 + d[5 ]);
      reg6  = v_load_64(score_matrix + 8 + d[6 ]);
      reg7  = v_load_64(score_matrix + 8 + d[7 ]);
      reg8  = v_load_64(score_matrix + 8 + d[8 ]);
      reg9  = v_load_64(score_matrix + 8 + d[9 ]);
      reg10 = v_load_64(score_matrix + 8 + d[10]);
      reg11 = v_load_64(score_matrix + 8 + d[11]);
      reg12 = v_load_64(score_matrix + 8 + d[12]);
      reg13 = v_load_64(score_matrix + 8 + d[13]);
      reg14 = v_load_64(score_matrix + 8 + d[14]);
      reg15 = v_load_64(score_matrix + 8 + d[15]);

      reg0  = v_merge_lo_8(reg0,  reg1);
      reg2  = v_merge_lo_8(reg2,  reg3);
      reg4  = v_merge_lo_8(reg4,  reg5);
      reg6  = v_merge_lo_8(reg6,  reg7);
      reg8  = v_merge_lo_8(reg8,  reg9);
      reg10 = v_merge_lo_8(reg10, reg11);
      reg12 = v_merge_lo_8(reg12, reg13);
      reg14 = v_merge_lo_8(reg14, reg15);

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile + 16 * j + 512 +   0, reg0);
      v_store(dprofile + 16 * j + 512 +  64, reg3);
      v_store(dprofile + 16 * j + 512 + 128, reg2);
      v_store(dprofile + 16 * j + 512 + 192, reg7);
      v_store(dprofile + 16 * j + 512 + 256, reg1);
      v_store(dprofile + 16 * j + 512 + 320, reg11);
      v_store(dprofile + 16 * j + 512 + 384, reg6);
      v_store(dprofile + 16 * j + 512 + 448, reg15);


      reg0  = v_load_64(score_matrix + 16 + d[0 ]);
      reg2  = v_load_64(score_matrix + 16 + d[2 ]);
      reg4  = v_load_64(score_matrix + 16 + d[4 ]);
      reg6  = v_load_64(score_matrix + 16 + d[6 ]);
      reg8  = v_load_64(score_matrix + 16 + d[8 ]);
      reg10 = v_load_64(score_matrix + 16 + d[10]);
      reg12 = v_load_64(score_matrix + 16 + d[12]);
      reg14 = v_load_64(score_matrix + 16 + d[14]);

      reg0  = v_merge_lo_8(reg0,  *CAST_VECTOR_p(score_matrix + 16 + d[ 1]));
      reg2  = v_merge_lo_8(reg2,  *CAST_VECTOR_p(score_matrix + 16 + d[ 3]));
      reg4  = v_merge_lo_8(reg4,  *CAST_VECTOR_p(score_matrix + 16 + d[ 5]));
      reg6  = v_merge_lo_8(reg6,  *CAST_VECTOR_p(score_matrix + 16 + d[ 7]));
      reg8  = v_merge_lo_8(reg8,  *CAST_VECTOR_p(score_matrix + 16 + d[ 9]));
      reg10 = v_merge_lo_8(reg10, *CAST_VECTOR_p(score_matrix + 16 + d[11]));
      reg12 = v_merge_lo_8(reg12, *CAST_VECTOR_p(score_matrix + 16 + d[13]));
      reg14 = v_merge_lo_8(reg14, *CAST_VECTOR_p(score_matrix + 16 + d[15]));

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile + 16 * j + 1024 +   0, reg0);
      v_store(dprofile + 16 * j + 1024 +  64, reg3);
      v_store(dprofile + 16 * j + 1024 + 128, reg2);
      v_store(dprofile + 16 * j + 1024 + 192, reg7);
      v_store(dprofile + 16 * j + 1024 + 256, reg1);
      v_store(dprofile + 16 * j + 1024 + 320, reg11);
      v_store(dprofile + 16 * j + 1024 + 384, reg6);
      v_store(dprofile + 16 * j + 1024 + 448, reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(score_matrix + 24 + d[ 0]);
      reg1  = v_load_64(score_matrix + 24 + d[ 1]);
      reg2  = v_load_64(score_matrix + 24 + d[ 2]);
      reg3  = v_load_64(score_matrix + 24 + d[ 3]);
      reg4  = v_load_64(score_matrix + 24 + d[ 4]);
      reg5  = v_load_64(score_matrix + 24 + d[ 5]);
      reg6  = v_load_64(score_matrix + 24 + d[ 6]);
      reg7  = v_load_64(score_matrix + 24 + d[ 7]);
      reg8  = v_load_64(score_matrix + 24 + d[ 8]);
      reg9  = v_load_64(score_matrix + 24 + d[ 9]);
      reg10 = v_load_64(score_matrix + 24 + d[10]);
      reg11 = v_load_64(score_matrix + 24 + d[11]);
      reg12 = v_load_64(score_matrix + 24 + d[12]);
      reg13 = v_load_64(score_matrix + 24 + d[13]);
      reg14 = v_load_64(score_matrix + 24 + d[14]);
      reg15 = v_load_64(score_matrix + 24 + d[15]);

      reg0  = v_merge_lo_8(reg0,  reg1);
      reg2  = v_merge_lo_8(reg2,  reg3);
      reg4  = v_merge_lo_8(reg4,  reg5);
      reg6  = v_merge_lo_8(reg6,  reg7);
      reg8  = v_merge_lo_8(reg8,  reg9);
      reg10 = v_merge_lo_8(reg10, reg11);
      reg12 = v_merge_lo_8(reg12, reg13);
      reg14 = v_merge_lo_8(reg14, reg15);

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile + 16 * j + 1536 +   0, reg0);  // refactoring: 1536 = 6 * 256 bits?
      v_store(dprofile + 16 * j + 1536 +  64, reg3);
      v_store(dprofile + 16 * j + 1536 + 128, reg2);
      v_store(dprofile + 16 * j + 1536 + 192, reg7);
      v_store(dprofile + 16 * j + 1536 + 256, reg1);
      v_store(dprofile + 16 * j + 1536 + 320, reg11);
      v_store(dprofile + 16 * j + 1536 + 384, reg6);
      v_store(dprofile + 16 * j + 1536 + 448, reg15);
    }
}

inline void onestep_8(VECTORTYPE & H,
                      VECTORTYPE & N,
                      VECTORTYPE & F,
                      VECTORTYPE V,
                      unsigned short * DIR,
                      VECTORTYPE & E,
                      VECTORTYPE QR,
                      VECTORTYPE R)
{
  VECTORTYPE W;

  H = v_add(H, V);
  W = H;
  H = v_min(H, F);
  *((DIR) + 0) = v_mask_eq(W, H);
  H = v_min(H, E);
  *((DIR) + 1) = v_mask_eq(H, E);
  N = H;
  H = v_add(H, QR);
  F = v_add(F, R);
  E = v_add(E, R);
  F = v_min(H, F);
  *((DIR) + 2) = v_mask_eq(H, F);
  E = v_min(H, E);
  *((DIR) + 3) = v_mask_eq(H, E);
}


void align_cells_regular_8(VECTORTYPE * Sm,
                           VECTORTYPE * hep,
                           VECTORTYPE ** qp,
                           VECTORTYPE * Qm,
                           VECTORTYPE * Rm,
                           uint64_t ql,
                           VECTORTYPE * F0,
                           uint64_t * dir_long,
                           VECTORTYPE * H0)
{
  VECTORTYPE Q;
  VECTORTYPE R;
  VECTORTYPE E;
  VECTORTYPE h0;
  VECTORTYPE h1;
  VECTORTYPE h2;
  VECTORTYPE h3;
  VECTORTYPE h4;
  VECTORTYPE h5;
  VECTORTYPE h6;
  VECTORTYPE h7;
  VECTORTYPE h8;
  VECTORTYPE f0;
  VECTORTYPE f1;
  VECTORTYPE f2;
  VECTORTYPE f3;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  Q = *Qm;
  R = *Rm;

  f0 = *F0;
  f1 = v_add(f0, R);
  f2 = v_add(f1, R);
  f3 = v_add(f2, R);

  h0 = *H0;
  h1 = v_sub(f0, Q);
  h2 = v_add(h1, R);
  h3 = v_add(h2, R);

  h5 = v_zero;
  h6 = v_zero;
  h7 = v_zero;
  h8 = v_zero;

  for(auto i = 0ULL; i < ql; i++)
    {
      VECTORTYPE *x {nullptr};
      static constexpr unsigned int j0 {0};
      static constexpr unsigned int j1 {j0 + 4};
      static constexpr unsigned int j2 {j1 + 4};
      static constexpr unsigned int j3 {j2 + 4};

      x = qp[i + 0];
      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];
      onestep_8(h0, h5, f0, x[0], dir + channels * i + j0, E, Q, R);
      onestep_8(h1, h6, f1, x[1], dir + channels * i + j1, E, Q, R);
      onestep_8(h2, h7, f2, x[2], dir + channels * i + j2, E, Q, R);
      onestep_8(h3, h8, f3, x[3], dir + channels * i + j3, E, Q, R);
      hep[2 * i + 0] = h8;
      hep[2 * i + 1] = E;
      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}


void align_cells_masked_8(VECTORTYPE * Sm,
                          VECTORTYPE * hep,
                          VECTORTYPE ** qp,
                          VECTORTYPE * Qm,
                          VECTORTYPE * Rm,
                          uint64_t ql,
                          VECTORTYPE * F0,
                          uint64_t * dir_long,
                          VECTORTYPE * H0,
                          VECTORTYPE * Mm,
                          VECTORTYPE * MQ,
                          VECTORTYPE * MR,
                          VECTORTYPE * MQ0)
{
  VECTORTYPE Q;
  VECTORTYPE R;
  VECTORTYPE E;
  VECTORTYPE h0;
  VECTORTYPE h1;
  VECTORTYPE h2;
  VECTORTYPE h3;
  VECTORTYPE h4;
  VECTORTYPE h5;
  VECTORTYPE h6;
  VECTORTYPE h7;
  VECTORTYPE h8;
  VECTORTYPE f0;
  VECTORTYPE f1;
  VECTORTYPE f2;
  VECTORTYPE f3;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  Q = *Qm;
  R = *Rm;

  f0 = *F0;
  f1 = v_add(f0, R);
  f2 = v_add(f1, R);
  f3 = v_add(f2, R);

  h0 = *H0;
  h1 = v_sub(f0, Q);
  h2 = v_add(h1, R);
  h3 = v_add(h2, R);

  h5 = v_zero;
  h6 = v_zero;
  h7 = v_zero;
  h8 = v_zero;

  for(auto i = 0ULL; i < ql; i++)
    {
      static constexpr unsigned int j0 {0};
      static constexpr unsigned int j1 {j0 + 4};
      static constexpr unsigned int j2 {j1 + 4};
      static constexpr unsigned int j3 {j2 + 4};

      VECTORTYPE * x {nullptr};

      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];
      x = qp[i + 0];

      /* mask h4 and E */
      h4 = v_sub(h4, *Mm);
      E  = v_sub(E,  *Mm);

      /* init h4 and E */
      h4 = v_add(h4, *MQ);
      E  = v_add(E,  *MQ);
      E  = v_add(E,  *MQ0);

      /* update MQ */
      *MQ = v_add(*MQ,  *MR);

      onestep_8(h0, h5, f0, x[0], dir + channels * i + j0, E, Q, R);
      onestep_8(h1, h6, f1, x[1], dir + channels * i + j1, E, Q, R);
      onestep_8(h2, h7, f2, x[2], dir + channels * i + j2, E, Q, R);
      onestep_8(h3, h8, f3, x[3], dir + channels * i + j3, E, Q, R);
      hep[2 * i + 0] = h8;
      hep[2 * i + 1] = E;

      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}


inline auto backtrack_8(char * qseq,
                        char * dseq,
                        uint64_t qlen,
                        uint64_t dlen,
                        uint64_t * dirbuffer,
                        uint64_t offset,
                        uint64_t dirbuffersize,
                        uint64_t channel,
                        uint64_t * alignmentlengthp,
                        const uint64_t longestdbsequence) -> uint64_t
{
  static constexpr unsigned int offset0 {0};
  static constexpr unsigned int offset1 {offset0 + 16};
  static constexpr unsigned int offset2 {offset1 + 16};
  static constexpr unsigned int offset3 {offset2 + 16};
  const uint64_t maskup      = 1ULL << (channel + offset0);
  const uint64_t maskleft    = 1ULL << (channel + offset1);
  const uint64_t maskextup   = 1ULL << (channel + offset2);
  const uint64_t maskextleft = 1ULL << (channel + offset3);

  auto i = static_cast<int64_t>(qlen) - 1;
  auto j = static_cast<int64_t>(dlen) - 1;
  uint64_t aligned {0};
  uint64_t matches {0};
  char operation {0};  // Insertion, Deletion or Match  // refactoring to enum class?

#undef SHOWALIGNMENT
#ifdef SHOWALIGNMENT
  printf("alignment, reversed: ");
#endif

  while ((i >= 0) and (j >= 0))
    {
      ++aligned;

      const uint64_t d
        = dirbuffer[(offset
                     + longestdbsequence * 4 * static_cast<uint64_t>(j / 4)
                     + static_cast<uint64_t>(4 * i + (j & 3U)))
                    % dirbuffersize];  // refactoring: how to rename that variable?

      if ((operation == 'I') and ((d & maskextleft) == 0U))
        {
          --j;
        }
      else if ((operation == 'D') and ((d & maskextup) == 0U))
        {
          --i;
        }
      else if ((d & maskleft) != 0U)
        {
          --j;
          operation = 'I';
        }
      else if ((d & maskup) == 0U)
        {
          --i;
          operation = 'D';
        }
      else
        {
          if (nt_extract(qseq, static_cast<uint64_t>(i)) ==
              nt_extract(dseq, static_cast<uint64_t>(j))) {
            ++matches;
          }
          --i;
          --j;
          operation = 'M';
        }

#ifdef SHOWALIGNMENT
      printf("%c", operation);
#endif
    }

  while (i >= 0)
    {
      ++aligned;
      --i;
#ifdef SHOWALIGNMENT
      printf("D");
#endif
    }

  while (j >= 0)
    {
      ++aligned;
      --j;
#ifdef SHOWALIGNMENT
      printf("I");
#endif
    }

#ifdef SHOWALIGNMENT
  printf("\n");
#endif

  * alignmentlengthp = aligned;
  return aligned - matches;
}


void search8(BYTE * * q_start,
             BYTE gap_open_penalty,
             BYTE gap_extend_penalty,
             BYTE * score_matrix,
             BYTE * dprofile,
             BYTE * hearray,
             uint64_t sequences,
             uint64_t * seqnos,
             uint64_t * scores,
             uint64_t * diffs,
             uint64_t * alignmentlengths,
             uint64_t qlen,
             uint64_t dirbuffersize,
             uint64_t * dirbuffer,
             const uint64_t longestdbsequence)
{
  VECTORTYPE Q;
  VECTORTYPE R;
  VECTORTYPE T;
  VECTORTYPE M;
  VECTORTYPE T0;
  VECTORTYPE MQ;
  VECTORTYPE MR;
  VECTORTYPE MQ0;
  VECTORTYPE *hep {nullptr};
  VECTORTYPE **qp {nullptr};

  std::array<uint64_t, channels> d_pos {{}};
  std::array<uint64_t, channels> d_offset {{}};
  std::array<char *, channels> d_address {{}};
  std::array<uint64_t, channels> d_length {{}};

  VECTORTYPE dseqalloc[cdepth];

  VECTORTYPE H0;
  VECTORTYPE F0;
  VECTORTYPE S[4];

  BYTE * dseq = reinterpret_cast<BYTE*>(& dseqalloc);

  std::array<int64_t, channels> seq_id {{}};
  uint64_t next_id {0};
  uint64_t done {0};

#ifdef __aarch64__
  const VECTORTYPE T0_init = { -1, 0, 0, 0, 0, 0, 0, 0,
                                0, 0, 0, 0, 0, 0, 0, 0 };
#elif defined __x86_64__
  const VECTORTYPE T0_init = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0,
                                          0, 0, 0, 0, 0, 0, 0, -1);
#elif defined __PPC__
  const VECTORTYPE T0_init = { (unsigned char)(-1), 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, 0 };
#endif

  T0 = T0_init;

  Q  = v_dup(static_cast<char>(gap_open_penalty + gap_extend_penalty));
  R  = v_dup(static_cast<char>(gap_extend_penalty));

  done = 0;

  hep = CAST_VECTOR_p(hearray);
  qp = reinterpret_cast<VECTORTYPE**>(q_start);

  // refactoring: initialize arrays
  for(auto c = 0U; c < channels; c++)
    {
      d_address[c] = nullptr;
      d_pos[c] = 0;
      d_length[c] = 0;
      seq_id[c] = -1;
    }

  F0 = v_zero;
  H0 = v_zero;

  bool easy {false};

  uint64_t * dir = dirbuffer;

  while(true)
    {
      if (easy)
        {
          // fill all channels

          for(auto c = 0U; c < channels; c++)
            {
              for(auto j = 0U; j < cdepth; j++)
                {
                  if (d_pos[c] < d_length[c]) {
                    dseq[channels * j + c]
                      = 1 + nt_extract(d_address[c], d_pos[c]++);
                  }
                  else {
                    dseq[channels * j + c] = 0;
                  }
                }
              if (d_pos[c] == d_length[c]) {
                easy = false;
              }
            }

#ifdef __x86_64__
          if (ssse3_present != 0)
            {
              dprofile_shuffle8(dprofile, score_matrix, dseq);
            }
          else
#endif
            {
              dprofile_fill8(dprofile, score_matrix, dseq);
            }

          align_cells_regular_8(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
        }
      else
        {
          // One or more sequences ended in the previous block
          // We have to switch over to a new sequence

          easy = true;

          M = v_zero;
          T = T0;
          for(auto c = 0U; c < channels; c++)
            {
              if (d_pos[c] < d_length[c])
                {
                  // this channel has more sequence

                  for(auto j = 0U; j < cdepth; j++)
                    {
                      if (d_pos[c] < d_length[c]) {
                        dseq[channels * j + c]
                          = 1 + nt_extract(d_address[c], d_pos[c]++);
                      }
                      else {
                        dseq[channels * j + c] = 0;
                      }
                    }
                  if (d_pos[c] == d_length[c]) {
                    easy = false;
                  }
                }
              else
                {
                  // sequence in channel c ended
                  // change of sequence

                  M = v_xor(M, T);

                  const int64_t cand_id = seq_id[c];

                  if (cand_id >= 0)
                    {
                      // save score

                      char * dbseq = d_address[c];
                      const uint64_t dbseqlen = d_length[c];
                      const uint64_t z = (dbseqlen + 3) % 4;
                      const uint64_t score
                        = (reinterpret_cast<BYTE*>(S))[z * channels + c];
                      scores[cand_id] = score;

                      uint64_t diff {0};

                      if (score < UINT8_MAX)
                        {
                          const uint64_t offset = d_offset[c];
                          diff = backtrack_8(query.seq, dbseq, qlen, dbseqlen,
                                             dirbuffer,
                                             offset,
                                             dirbuffersize, c,
                                             alignmentlengths + cand_id,
                                             longestdbsequence);
                        }
                      else
                        {
                          diff = UINT8_MAX;
                        }

                      diffs[cand_id] = diff;

                      ++done;
                    }

                  if (next_id < sequences)
                    {
                      // get next sequence
                      seq_id[c] = static_cast<int64_t>(next_id);
                      const uint64_t seqno = seqnos[next_id];
                      char * address {nullptr};
                      unsigned int length {0};

                      db_getsequenceandlength(seqno, & address, & length);

                      d_address[c] = address;
                      d_length[c] = length;

                      d_pos[c] = 0;
                      d_offset[c] = static_cast<uint64_t>(dir - dirbuffer);
                      ++next_id;

                      (reinterpret_cast<BYTE*>(&H0))[c] = 0;
                      (reinterpret_cast<BYTE*>(&F0))[c] = static_cast<BYTE>(2U * gap_open_penalty + 2U * gap_extend_penalty);

                      // fill channel
                      for(auto j = 0U; j < cdepth; j++)
                        {
                          if (d_pos[c] < d_length[c]) {
                            dseq[channels * j + c] = 1 + nt_extract(d_address[c], d_pos[c]++);
                          }
                          else {
                            dseq[channels * j + c] = 0;
                          }
                        }
                      if (d_pos[c] == d_length[c]) {
                        easy = false;
                      }
                    }
                  else
                    {
                      // no more sequences, empty channel
                      seq_id[c] = -1;
                      d_address[c] = nullptr;
                      d_pos[c] = 0;
                      d_length[c] = 0;
                      for(auto j = 0U; j < cdepth; j++) {
                        dseq[channels * j + c] = 0;
                      }
                    }
                }

              T = v_shift_left(T);
            }

          if (done == sequences) {
            break;
          }

#ifdef __x86_64__
          if (ssse3_present != 0)
            {
              dprofile_shuffle8(dprofile, score_matrix, dseq);
            }
          else
#endif
            {
              dprofile_fill8(dprofile, score_matrix, dseq);
            }

          MQ = v_and(M, Q);
          MR = v_and(M, R);
          MQ0 = MQ;

          align_cells_masked_8(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
        }

      F0 = v_add(F0, R);
      F0 = v_add(F0, R);
      F0 = v_add(F0, R);
      H0 = v_sub(F0, Q);
      F0 = v_add(F0, R);

      dir += 4 * longestdbsequence;
      if (dir >= dirbuffer + dirbuffersize) {
        dir -= dirbuffersize;
      }
    }
}
