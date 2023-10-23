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

#ifdef __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#endif

#include "ssse3.h"
#include "sse41.h"

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


constexpr unsigned int channels {8};
constexpr unsigned int cdepth {4};
using BYTE = unsigned char;
using WORD = unsigned short;


/* uses 8 unsigned 16-bit values */

#ifdef __aarch64__

using VECTORTYPE = int16x8_t;

#define CAST_VECTOR_p(x) (reinterpret_cast<VECTORTYPE *>(x))

const uint16x8_t neon_mask =
  {0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000};

#define v_load(a) vld1q_s16((const int16_t *)(a))
#define v_store(a, b) vst1q_s16((int16_t *)(a), (b))
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
#define v_min(a, b) vminq_u16((a), (b))
#define v_add(a, b) vqaddq_u16((a), (b))
#define v_sub(a, b) vqsubq_u16((a), (b))
#define v_dup(a) vdupq_n_s16(a)
#define v_zero v_dup(0)
#define v_and(a, b) vandq_u16((a), (b))
#define v_xor(a, b) veorq_u16((a), (b))
#define v_shift_left(a) vextq_u16((v_zero), (a), 7)
#define v_mask_eq(a, b) vaddvq_u16(vandq_u16((vceqq_s16((a), (b))), neon_mask))

#elif defined __x86_64__

using VECTORTYPE = __m128i;

#define CAST_VECTOR_p(x) (reinterpret_cast<VECTORTYPE *>(x))

#define v_load(a) _mm_load_si128(CAST_VECTOR_p(a))
#define v_store(a, b) _mm_store_si128(CAST_VECTOR_p(a), (b))
#define v_merge_lo_16(a, b) _mm_unpacklo_epi16((a),(b))
#define v_merge_hi_16(a, b) _mm_unpackhi_epi16((a),(b))
#define v_merge_lo_32(a, b) _mm_unpacklo_epi32((a),(b))
#define v_merge_hi_32(a, b) _mm_unpackhi_epi32((a),(b))
#define v_merge_lo_64(a, b) _mm_unpacklo_epi64((a),(b))
#define v_merge_hi_64(a, b) _mm_unpackhi_epi64((a),(b))
#define v_min(a, b) _mm_subs_epu16((a), _mm_subs_epu16((a), (b)))
#define v_add(a, b) _mm_adds_epu16((a), (b))
#define v_sub(a, b) _mm_subs_epu16((a), (b))
#define v_dup(a) _mm_set1_epi16(a)
#define v_zero v_dup(0)
#define v_and(a, b) _mm_and_si128((a), (b))
#define v_xor(a, b) _mm_xor_si128((a), (b))
#define v_shift_left(a) _mm_slli_si128((a), 2)
#define v_mask_eq(a, b) static_cast<unsigned short> \
  (_mm_movemask_epi8(_mm_cmpeq_epi16((a), (b))))

#elif defined __PPC__

using VECTORTYPE = vector unsigned short;

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

#define v_load(a) *CAST_VECTOR_p(a)
#define v_store(a, b) vec_st((VECTORTYPE)(b), 0, CAST_VECTOR_p(a))
#define v_merge_lo_16(a, b) vec_mergeh((VECTORTYPE)(a), (VECTORTYPE)(b))
#define v_merge_hi_16(a, b) vec_mergel((VECTORTYPE)(a), (VECTORTYPE)(b))
#define v_merge_lo_32(a, b) (VECTORTYPE) vec_mergeh((vector int)(a),    \
                                                    (vector int)(b))
#define v_merge_hi_32(a, b) (VECTORTYPE) vec_mergel((vector int)(a),    \
                                                    (vector int)(b))
#define v_merge_lo_64(a, b) (VECTORTYPE) vec_perm((vector long long)(a), \
                                                  (vector long long)(b), \
                                                  perm_merge_long_low)
#define v_merge_hi_64(a, b) (VECTORTYPE) vec_perm((vector long long)(a), \
                                                  (vector long long)(b), \
                                                  perm_merge_long_high)
#define v_min(a, b) vec_min((a), (b))
#define v_add(a, b) vec_adds((a), (b))
#define v_sub(a, b) vec_subs((a), (b))
#define v_dup(a) vec_splats((unsigned short)(a));
#define v_zero vec_splat_u16(0)
#define v_and(a, b) vec_and((a), (b))
#define v_xor(a, b) vec_xor((a), (b))
#define v_shift_left(a) vec_sld((a), v_zero, 2)
#define v_mask_eq(a, b) ((vector unsigned short) \
  vec_vbpermq((vector unsigned char) vec_cmpeq((a), (b)), perm_bits))[4]

#else

#error Unknown Architecture

#endif


inline void dprofile_fill16(WORD * dprofile_word,
                            WORD * score_matrix_word,
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
  VECTORTYPE reg16;
  VECTORTYPE reg17;
  VECTORTYPE reg18;
  VECTORTYPE reg19;
  VECTORTYPE reg20;
  VECTORTYPE reg21;
  VECTORTYPE reg22;
  VECTORTYPE reg23;
  VECTORTYPE reg24;
  VECTORTYPE reg25;
  VECTORTYPE reg26;
  VECTORTYPE reg27;
  VECTORTYPE reg28;
  VECTORTYPE reg29;
  VECTORTYPE reg30;
  VECTORTYPE reg31;

  for(auto j = 0U; j < cdepth; j++)
    {
      std::array<unsigned int, channels> d {{}};
      for(auto z = 0U; z < channels; z++) {
        d[z] = (static_cast<unsigned int>(dseq[j * channels + z])) << multiplier;
      }

      for(auto i = 0U; i < channels; i += channels)
        {
          reg0  = v_load(score_matrix_word + d[0] + i);
          reg1  = v_load(score_matrix_word + d[1] + i);
          reg2  = v_load(score_matrix_word + d[2] + i);
          reg3  = v_load(score_matrix_word + d[3] + i);
          reg4  = v_load(score_matrix_word + d[4] + i);
          reg5  = v_load(score_matrix_word + d[5] + i);
          reg6  = v_load(score_matrix_word + d[6] + i);
          reg7  = v_load(score_matrix_word + d[7] + i);

          reg8  = v_merge_lo_16(reg0,  reg1);
          reg9  = v_merge_hi_16(reg0,  reg1);
          reg10 = v_merge_lo_16(reg2,  reg3);
          reg11 = v_merge_hi_16(reg2,  reg3);
          reg12 = v_merge_lo_16(reg4,  reg5);
          reg13 = v_merge_hi_16(reg4,  reg5);
          reg14 = v_merge_lo_16(reg6,  reg7);
          reg15 = v_merge_hi_16(reg6,  reg7);

          reg16 = v_merge_lo_32(reg8,  reg10);
          reg17 = v_merge_hi_32(reg8,  reg10);
          reg18 = v_merge_lo_32(reg12, reg14);
          reg19 = v_merge_hi_32(reg12, reg14);
          reg20 = v_merge_lo_32(reg9,  reg11);
          reg21 = v_merge_hi_32(reg9,  reg11);
          reg22 = v_merge_lo_32(reg13, reg15);
          reg23 = v_merge_hi_32(reg13, reg15);

          reg24 = v_merge_lo_64(reg16, reg18);
          reg25 = v_merge_hi_64(reg16, reg18);
          reg26 = v_merge_lo_64(reg17, reg19);
          reg27 = v_merge_hi_64(reg17, reg19);
          reg28 = v_merge_lo_64(reg20, reg22);
          reg29 = v_merge_hi_64(reg20, reg22);
          reg30 = v_merge_lo_64(reg21, reg23);
          reg31 = v_merge_hi_64(reg21, reg23);

          v_store(dprofile_word + cdepth * channels * (i + 0) + channels * j, reg24);
          v_store(dprofile_word + cdepth * channels * (i + 1) + channels * j, reg25);
          v_store(dprofile_word + cdepth * channels * (i + 2) + channels * j, reg26);
          v_store(dprofile_word + cdepth * channels * (i + 3) + channels * j, reg27);
          v_store(dprofile_word + cdepth * channels * (i + 4) + channels * j, reg28);
          v_store(dprofile_word + cdepth * channels * (i + 5) + channels * j, reg29);
          v_store(dprofile_word + cdepth * channels * (i + 6) + channels * j, reg30);
          v_store(dprofile_word + cdepth * channels * (i + 7) + channels * j, reg31);
        }
    }
}

inline void onestep_16(VECTORTYPE & H,
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
  *(DIR+0) = v_mask_eq(W, H);
  H = v_min(H, E);
  *(DIR+1) = v_mask_eq(H, E);
  N = H;
  H = v_add(H, QR);
  F = v_add(F, R);
  E = v_add(E, R);
  F = v_min(H, F);
  *(DIR+2) = v_mask_eq(H, F);
  E = v_min(H, E);
  *(DIR+3) = v_mask_eq(H, E);
}

void align_cells_regular_16(VECTORTYPE * Sm,
                            VECTORTYPE * hep,
                            VECTORTYPE ** qp,
                            VECTORTYPE * Qm,
                            VECTORTYPE * Rm,
                            uint64_t ql,
                            VECTORTYPE * F0,
                            uint64_t * dir_long,
                            VECTORTYPE * H0)
{
  static constexpr unsigned int step {16};
  static constexpr unsigned int offset0 {0};
  static constexpr unsigned int offset1 {offset0 + 4};
  static constexpr unsigned int offset2 {offset1 + 4};
  static constexpr unsigned int offset3 {offset2 + 4};

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

      x = qp[i + 0];
      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];
      onestep_16(h0, h5, f0, x[0], dir + step * i + offset0, E, Q, R);
      onestep_16(h1, h6, f1, x[1], dir + step * i + offset1, E, Q, R);
      onestep_16(h2, h7, f2, x[2], dir + step * i + offset2, E, Q, R);
      onestep_16(h3, h8, f3, x[3], dir + step * i + offset3, E, Q, R);
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

void align_cells_masked_16(VECTORTYPE * Sm,
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
  static constexpr unsigned int step {16};
  static constexpr unsigned int offset0 {0};
  static constexpr unsigned int offset1 {offset0 + 4};
  static constexpr unsigned int offset2 {offset1 + 4};
  static constexpr unsigned int offset3 {offset2 + 4};

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

      onestep_16(h0, h5, f0, x[0], dir + step * i + offset0, E, Q, R);
      onestep_16(h1, h6, f1, x[1], dir + step * i + offset1, E, Q, R);
      onestep_16(h2, h7, f2, x[2], dir + step * i + offset2, E, Q, R);
      onestep_16(h3, h8, f3, x[3], dir + step * i + offset3, E, Q, R);
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

auto backtrack_16(char * qseq,
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
  const uint64_t maskup      = 3ULL << (2 * channel + offset0);
  const uint64_t maskleft    = 3ULL << (2 * channel + offset1);
  const uint64_t maskextup   = 3ULL << (2 * channel + offset2);
  const uint64_t maskextleft = 3ULL << (2 * channel + offset3);

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

void search16(WORD * * q_start,
              WORD gap_open_penalty,
              WORD gap_extend_penalty,
              WORD * score_matrix,
              WORD * dprofile,
              WORD * hearray,
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
  const VECTORTYPE T0_init = { -1, 0, 0, 0, 0, 0, 0, 0 };
#elif defined __x86_64__
  const VECTORTYPE T0_init = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, -1);
#elif defined __PPC__
  const VECTORTYPE T0_init = { (unsigned short)(-1), 0, 0, 0, 0, 0, 0, 0 };
#endif

  T0 = T0_init;

  Q  = v_dup(static_cast<short>(gap_open_penalty + gap_extend_penalty));
  R  = v_dup(static_cast<short>(gap_extend_penalty));

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
              dprofile_shuffle16(dprofile, score_matrix, dseq);
            }
          else
#endif
            {
              dprofile_fill16(dprofile, score_matrix, dseq);
            }

#ifdef __x86_64__
          if (sse41_present != 0)
            {
              align_cells_regular_16_sse41(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
            }
          else
#endif
            {
              align_cells_regular_16(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
            }
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
                        dseq[channels*j+c] = 0;
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

                      char * dbseq = reinterpret_cast<char*>(d_address[c]);
                      const uint64_t dbseqlen = d_length[c];
                      const uint64_t z = (dbseqlen + 3) % 4;
                      const uint64_t score
                        = (reinterpret_cast<WORD*>(S))[z * channels + c];
                      scores[cand_id] = score;

                      uint64_t diff {0};

                      if (score < UINT16_MAX)
                        {
                          const uint64_t offset = d_offset[c];
                          diff = backtrack_16(query.seq, dbseq, qlen, dbseqlen,
                                              dirbuffer,
                                              offset,
                                              dirbuffersize, c,
                                              alignmentlengths + cand_id,
                                              longestdbsequence);
                        }
                      else
                        {
                          diff = UINT16_MAX;
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

                      (reinterpret_cast<WORD*>(&H0))[c] = 0;
                      (reinterpret_cast<WORD*>(&F0))[c] = static_cast<WORD>(2U * gap_open_penalty + 2U * gap_extend_penalty);

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
              dprofile_shuffle16(dprofile, score_matrix, dseq);
            }
          else
#endif
            {
              dprofile_fill16(dprofile, score_matrix, dseq);
            }

          MQ = v_and(M, Q);
          MR = v_and(M, R);
          MQ0 = MQ;

#ifdef __x86_64__
          if (sse41_present != 0)
            {
              align_cells_masked_16_sse41(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
            }
          else
#endif
            {
              align_cells_masked_16(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
            }
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
