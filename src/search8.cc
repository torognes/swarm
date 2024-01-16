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

#include "swarm.h"
#include "db.h"
#include "utils/backtrack.h"
#include <array>
#include <cstdint>  // int64_t, uint64_t, uint8_t
#include <limits>


// refactoring: C++26 std::simd
#ifdef __aarch64__

#include <arm_neon.h>
#include "utils/intrinsics_to_functions_aarch64.h"
using VECTORTYPE = uint8x16_t;

#elif defined __x86_64__

#ifdef __SSE2__

#include <emmintrin.h>  // SSE2 intrinsics
#include "utils/intrinsics_to_functions_x86_64.h"
using VECTORTYPE = __m128i;

#endif

#include "ssse3.h"

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__

#include <altivec.h>
#include "utils/intrinsics_to_functions_ppc.h"
using VECTORTYPE = vector unsigned char;

#else

#error Big endian ppc64 CPUs not supported

#endif

#else

#error Unknown architecture

#endif


constexpr unsigned int channels {16};
constexpr unsigned int cdepth {4};
constexpr uint8_t n_bits {8};
using BYTE = unsigned char;

// backtrack.h: template specialization (8 bits)
template <>
auto compute_mask<n_bits>(uint64_t const channel,
                     unsigned int const offset) -> uint64_t {
  return (1ULL << (channel + offset));
}

// refactoring: objdump shows this function is not inlined
inline void dprofile_fill8(BYTE * dprofile,
                           BYTE * score_matrix,
                           BYTE const * dseq)
{
  static constexpr auto multiplier = 5U;
  static constexpr auto n_lanes = 16ULL;  // refactoring: same as channels?

  static constexpr auto pos0 = 0U;
  static constexpr auto pos1 = pos0 + 1;
  static constexpr auto pos2 = pos1 + 1;
  static constexpr auto pos3 = pos2 + 1;
  static constexpr auto pos4 = pos3 + 1;
  static constexpr auto pos5 = pos4 + 1;
  static constexpr auto pos6 = pos5 + 1;
  static constexpr auto pos7 = pos6 + 1;
  static constexpr auto pos8 = pos7 + 1;
  static constexpr auto pos9 = pos8 + 1;
  static constexpr auto pos10 = pos9 + 1;
  static constexpr auto pos11 = pos10 + 1;
  static constexpr auto pos12 = pos11 + 1;
  static constexpr auto pos13 = pos12 + 1;
  static constexpr auto pos14 = pos13 + 1;
  static constexpr auto pos15 = pos14 + 1;

  static constexpr auto line0 = 64U * 0;  // as in 'cache line': 64 bytes
  static constexpr auto line1 = 64U * 1;
  static constexpr auto line2 = 64U * 2;
  static constexpr auto line3 = 64U * 3;
  static constexpr auto line4 = 64U * 4;
  static constexpr auto line5 = 64U * 5;
  static constexpr auto line6 = 64U * 6;
  static constexpr auto line7 = 64U * 7;
  static constexpr auto line8 = 64U * 8;
  static constexpr auto line16 = 64U * 16;  // 1,024
  static constexpr auto line24 = 64U * 24;  // 1,536

  static constexpr auto offset8 = 8U;
  static constexpr auto offset16 = 16U;
  static constexpr auto offset24 = 24U;

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

      reg0  = v_load_64(score_matrix + d[pos0]);
      reg2  = v_load_64(score_matrix + d[pos2]);
      reg4  = v_load_64(score_matrix + d[pos4]);
      reg6  = v_load_64(score_matrix + d[pos6]);
      reg8  = v_load_64(score_matrix + d[pos8]);
      reg10 = v_load_64(score_matrix + d[pos10]);
      reg12 = v_load_64(score_matrix + d[pos12]);
      reg14 = v_load_64(score_matrix + d[pos14]);

      reg0  = v_merge_lo_8(reg0,  *cast_vector8(score_matrix + d[pos1]));
      reg2  = v_merge_lo_8(reg2,  *cast_vector8(score_matrix + d[pos3]));
      reg4  = v_merge_lo_8(reg4,  *cast_vector8(score_matrix + d[pos5]));
      reg6  = v_merge_lo_8(reg6,  *cast_vector8(score_matrix + d[pos7]));
      reg8  = v_merge_lo_8(reg8,  *cast_vector8(score_matrix + d[pos9]));
      reg10 = v_merge_lo_8(reg10, *cast_vector8(score_matrix + d[pos11]));
      reg12 = v_merge_lo_8(reg12, *cast_vector8(score_matrix + d[pos13]));
      reg14 = v_merge_lo_8(reg14, *cast_vector8(score_matrix + d[pos15]));

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

      v_store8(cast_vector8(dprofile + n_lanes * j + line0), reg0);
      v_store8(cast_vector8(dprofile + n_lanes * j + line1), reg3);
      v_store8(cast_vector8(dprofile + n_lanes * j + line2), reg2);
      v_store8(cast_vector8(dprofile + n_lanes * j + line3), reg7);
      v_store8(cast_vector8(dprofile + n_lanes * j + line4), reg1);
      v_store8(cast_vector8(dprofile + n_lanes * j + line5), reg11);
      v_store8(cast_vector8(dprofile + n_lanes * j + line6), reg6);
      v_store8(cast_vector8(dprofile + n_lanes * j + line7), reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(score_matrix + offset8 + d[pos0]);
      reg1  = v_load_64(score_matrix + offset8 + d[pos1]);
      reg2  = v_load_64(score_matrix + offset8 + d[pos2]);
      reg3  = v_load_64(score_matrix + offset8 + d[pos3]);
      reg4  = v_load_64(score_matrix + offset8 + d[pos4]);
      reg5  = v_load_64(score_matrix + offset8 + d[pos5]);
      reg6  = v_load_64(score_matrix + offset8 + d[pos6]);
      reg7  = v_load_64(score_matrix + offset8 + d[pos7]);
      reg8  = v_load_64(score_matrix + offset8 + d[pos8]);
      reg9  = v_load_64(score_matrix + offset8 + d[pos9]);
      reg10 = v_load_64(score_matrix + offset8 + d[pos10]);
      reg11 = v_load_64(score_matrix + offset8 + d[pos11]);
      reg12 = v_load_64(score_matrix + offset8 + d[pos12]);
      reg13 = v_load_64(score_matrix + offset8 + d[pos13]);
      reg14 = v_load_64(score_matrix + offset8 + d[pos14]);
      reg15 = v_load_64(score_matrix + offset8 + d[pos15]);

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

      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line0), reg0);
      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line1), reg3);
      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line2), reg2);
      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line3), reg7);
      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line4), reg1);
      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line5), reg11);
      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line6), reg6);
      v_store8(cast_vector8(dprofile + n_lanes * j + line8 + line7), reg15);


      reg0  = v_load_64(score_matrix + offset16 + d[pos0]);
      reg2  = v_load_64(score_matrix + offset16 + d[pos2]);
      reg4  = v_load_64(score_matrix + offset16 + d[pos4]);
      reg6  = v_load_64(score_matrix + offset16 + d[pos6]);
      reg8  = v_load_64(score_matrix + offset16 + d[pos8]);
      reg10 = v_load_64(score_matrix + offset16 + d[pos10]);
      reg12 = v_load_64(score_matrix + offset16 + d[pos12]);
      reg14 = v_load_64(score_matrix + offset16 + d[pos14]);

      reg0  = v_merge_lo_8(reg0,  *cast_vector8(score_matrix + offset16 + d[pos1]));
      reg2  = v_merge_lo_8(reg2,  *cast_vector8(score_matrix + offset16 + d[pos3]));
      reg4  = v_merge_lo_8(reg4,  *cast_vector8(score_matrix + offset16 + d[pos5]));
      reg6  = v_merge_lo_8(reg6,  *cast_vector8(score_matrix + offset16 + d[pos7]));
      reg8  = v_merge_lo_8(reg8,  *cast_vector8(score_matrix + offset16 + d[pos9]));
      reg10 = v_merge_lo_8(reg10, *cast_vector8(score_matrix + offset16 + d[pos11]));
      reg12 = v_merge_lo_8(reg12, *cast_vector8(score_matrix + offset16 + d[pos13]));
      reg14 = v_merge_lo_8(reg14, *cast_vector8(score_matrix + offset16 + d[pos15]));

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

      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line0), reg0);
      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line1), reg3);
      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line2), reg2);
      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line3), reg7);
      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line4), reg1);
      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line5), reg11);
      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line6), reg6);
      v_store8(cast_vector8(dprofile + n_lanes * j + line16 + line7), reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(score_matrix + offset24 + d[pos0]);
      reg1  = v_load_64(score_matrix + offset24 + d[pos1]);
      reg2  = v_load_64(score_matrix + offset24 + d[pos2]);
      reg3  = v_load_64(score_matrix + offset24 + d[pos3]);
      reg4  = v_load_64(score_matrix + offset24 + d[pos4]);
      reg5  = v_load_64(score_matrix + offset24 + d[pos5]);
      reg6  = v_load_64(score_matrix + offset24 + d[pos6]);
      reg7  = v_load_64(score_matrix + offset24 + d[pos7]);
      reg8  = v_load_64(score_matrix + offset24 + d[pos8]);
      reg9  = v_load_64(score_matrix + offset24 + d[pos9]);
      reg10 = v_load_64(score_matrix + offset24 + d[pos10]);
      reg11 = v_load_64(score_matrix + offset24 + d[pos11]);
      reg12 = v_load_64(score_matrix + offset24 + d[pos12]);
      reg13 = v_load_64(score_matrix + offset24 + d[pos13]);
      reg14 = v_load_64(score_matrix + offset24 + d[pos14]);
      reg15 = v_load_64(score_matrix + offset24 + d[pos15]);

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

      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line0), reg0);
      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line1), reg3);
      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line2), reg2);
      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line3), reg7);
      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line4), reg1);
      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line5), reg11);
      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line6), reg6);
      v_store8(cast_vector8(dprofile + n_lanes * j + line24 + line7), reg15);
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
  H = v_add8(H, V);
  const auto W = H;
  H = v_min8(H, F);
  *((DIR) + 0) = v_mask_eq8(W, H);
  H = v_min8(H, E);
  *((DIR) + 1) = v_mask_eq8(H, E);
  N = H;
  H = v_add8(H, QR);
  F = v_add8(F, R);
  E = v_add8(E, R);
  F = v_min8(H, F);
  *((DIR) + 2) = v_mask_eq8(H, F);
  E = v_min8(H, E);
  *((DIR) + 3) = v_mask_eq8(H, E);
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
  static constexpr auto offset0 = 0U;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  const auto Q = *Qm;
  const auto R = *Rm;

  auto f0 = *F0;
  auto f1 = v_add8(f0, R);
  auto f2 = v_add8(f1, R);
  auto f3 = v_add8(f2, R);

  auto h0 = *H0;
  auto h1 = v_sub8(f0, Q);
  auto h2 = v_add8(h1, R);
  auto h3 = v_add8(h2, R);

  auto h5 = v_zero8();
  auto h6 = v_zero8();
  auto h7 = v_zero8();
  auto h8 = v_zero8();

  for(auto i = 0ULL; i < ql; i++)
    {
      VECTORTYPE *x = qp[i + 0];
      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];
      onestep_8(h0, h5, f0, x[0], dir + channels * i + offset0, E, Q, R);
      onestep_8(h1, h6, f1, x[1], dir + channels * i + offset1, E, Q, R);
      onestep_8(h2, h7, f2, x[2], dir + channels * i + offset2, E, Q, R);
      onestep_8(h3, h8, f3, x[3], dir + channels * i + offset3, E, Q, R);
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
  static constexpr auto offset0 = 0U;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * dir = reinterpret_cast<unsigned short *>(dir_long);

  const auto Q = *Qm;
  const auto R = *Rm;

  auto f0 = *F0;
  auto f1 = v_add8(f0, R);
  auto f2 = v_add8(f1, R);
  auto f3 = v_add8(f2, R);

  auto h0 = *H0;
  auto h1 = v_sub8(f0, Q);
  auto h2 = v_add8(h1, R);
  auto h3 = v_add8(h2, R);

  auto h5 = v_zero8();
  auto h6 = v_zero8();
  auto h7 = v_zero8();
  auto h8 = v_zero8();

  for(auto i = 0ULL; i < ql; i++)
    {
      VECTORTYPE * x = qp[i + 0];

      h4 = hep[2 * i + 0];
      E  = hep[2 * i + 1];

      /* mask h4 and E */
      h4 = v_sub8(h4, *Mm);
      E  = v_sub8(E,  *Mm);

      /* init h4 and E */
      h4 = v_add8(h4, *MQ);
      E  = v_add8(E,  *MQ);
      E  = v_add8(E,  *MQ0);

      /* update MQ */
      *MQ = v_add8(*MQ,  *MR);

      onestep_8(h0, h5, f0, x[0], dir + channels * i + offset0, E, Q, R);
      onestep_8(h1, h6, f1, x[1], dir + channels * i + offset1, E, Q, R);
      onestep_8(h2, h7, f2, x[2], dir + channels * i + offset2, E, Q, R);
      onestep_8(h3, h8, f3, x[3], dir + channels * i + offset3, E, Q, R);
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


auto search8(BYTE * * q_start,
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
             const uint64_t longestdbsequence) -> void
{
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
  VECTORTYPE T;
  VECTORTYPE M;
  VECTORTYPE MQ;
  VECTORTYPE MR;
  VECTORTYPE MQ0;

  // by default, std::array is value-initialized (set to 0 for int,
  // nullptr for pointers, etc)
  std::array<uint64_t, channels> d_pos {{}};
  std::array<uint64_t, channels> d_offset {{}};
  std::array<char *, channels> d_address {{}};
  std::array<uint64_t, channels> d_length {{}};
  std::array<int64_t, channels> seq_id {{}};
  seq_id.fill(-1);

  // refactoring fail: std::array -> warning: ignoring attributes on
  // template argument ‘VECTORTYPE’ {aka ‘__m128i’}
  VECTORTYPE S[4];

  // make an array of size VECTORTYPE * channels, but interpret as
  // an array of BYTES
  std::array<BYTE, channels * sizeof(VECTORTYPE) / sizeof(BYTE)> dseq {{}};

  uint64_t next_id {0};
  uint64_t done {0};

#ifdef __aarch64__
  const VECTORTYPE T0 = { uint8_max, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0 };
#elif defined __x86_64__
  const auto T0 = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0,
                               0, 0, 0, 0, 0, 0, 0, -1);
#elif defined __PPC__
  static constexpr auto uchar_max = std::numeric_limits<unsigned char>::max();
  const VECTORTYPE T0 = { uchar_max, 0, 0, 0, 0, 0, 0, 0,
                                  0, 0, 0, 0, 0, 0, 0, 0 };
#endif

  auto Q = v_dup8(static_cast<char>(gap_open_penalty + gap_extend_penalty));
  auto R = v_dup8(static_cast<char>(gap_extend_penalty));

  done = 0;

  // refactoring: can't remove reinterpret_cast, cast_vector8() is a nullop in Aarch64
  auto *hep = reinterpret_cast<VECTORTYPE*>(hearray);
  auto **qp = reinterpret_cast<VECTORTYPE**>(q_start);

  auto F0 = v_zero8();
  auto H0 = v_zero8();

  bool easy {false};

  uint64_t * dir = dirbuffer;

  while(true)
    {
      if (easy)
        {
          // fill all channels

          for(auto channel = 0U; channel < channels; channel++)
            {
              for(auto j = 0U; j < cdepth; j++)
                {
                  if (d_pos[channel] < d_length[channel]) {
                    dseq[channels * j + channel]
                      = 1 + nt_extract(d_address[channel], d_pos[channel]++);
                  }
                  else {
                    dseq[channels * j + channel] = 0;
                  }
                }
              if (d_pos[channel] == d_length[channel]) {
                easy = false;
              }
            }

#ifdef __x86_64__
#ifdef __SSE3__
          if (ssse3_present != 0)
            {
              dprofile_shuffle8(dprofile, score_matrix, dseq.data());
            }
          else
#endif
#endif
            {
              dprofile_fill8(dprofile, score_matrix, dseq.data());
            }

          align_cells_regular_8(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
        }
      else
        {
          // One or more sequences ended in the previous block
          // We have to switch over to a new sequence

          easy = true;

          M = v_zero8();
          T = T0;
          for(auto channel = 0U; channel < channels; channel++)
            {
              if (d_pos[channel] < d_length[channel])
                {
                  // this channel has more sequence

                  for(auto j = 0U; j < cdepth; j++)
                    {
                      if (d_pos[channel] < d_length[channel]) {
                        dseq[channels * j + channel]
                          = 1 + nt_extract(d_address[channel], d_pos[channel]++);
                      }
                      else {
                        dseq[channels * j + channel] = 0;
                      }
                    }
                  if (d_pos[channel] == d_length[channel]) {
                    easy = false;
                  }
                }
              else
                {
                  // sequence in channel ended,
                  // change of sequence

                  M = v_xor8(M, T);

                  const int64_t cand_id = seq_id[channel];

                  if (cand_id >= 0)
                    {
                      // save score

                      char * dbseq = d_address[channel];
                      const uint64_t dbseqlen = d_length[channel];
                      const uint64_t z = (dbseqlen + 3) % 4;
                      const uint64_t score
                        = (reinterpret_cast<BYTE*>(S))[z * channels + channel];
                      scores[cand_id] = score;

                      uint64_t diff {0};

                      if (score < uint8_max)
                        {
                          const uint64_t offset = d_offset[channel];
                          diff = backtrack<n_bits>(query.seq, dbseq, qlen, dbseqlen,
                                                   dirbuffer,
                                                   offset,
                                                   dirbuffersize, channel,
                                                   alignmentlengths + cand_id,
                                                   longestdbsequence);
                        }
                      else
                        {
                          diff = uint8_max;
                        }

                      diffs[cand_id] = diff;

                      ++done;
                    }

                  if (next_id < sequences)
                    {
                      // get next sequence
                      seq_id[channel] = static_cast<int64_t>(next_id);
                      const uint64_t seqno = seqnos[next_id];
                      char * address {nullptr};
                      unsigned int length {0};

                      db_getsequenceandlength(seqno, & address, length);

                      d_address[channel] = address;
                      d_length[channel] = length;

                      d_pos[channel] = 0;
                      d_offset[channel] = static_cast<uint64_t>(dir - dirbuffer);
                      ++next_id;

                      (reinterpret_cast<BYTE*>(&H0))[channel] = 0;
                      (reinterpret_cast<BYTE*>(&F0))[channel] = static_cast<BYTE>(2U * gap_open_penalty + 2U * gap_extend_penalty);

                      // fill channel
                      for(auto j = 0U; j < cdepth; j++)
                        {
                          if (d_pos[channel] < d_length[channel]) {
                            dseq[channels * j + channel] = 1 + nt_extract(d_address[channel], d_pos[channel]++);
                          }
                          else {
                            dseq[channels * j + channel] = 0;
                          }
                        }
                      if (d_pos[channel] == d_length[channel]) {
                        easy = false;
                      }
                    }
                  else
                    {
                      // no more sequences, empty channel
                      seq_id[channel] = -1;
                      d_address[channel] = nullptr;
                      d_pos[channel] = 0;
                      d_length[channel] = 0;
                      for(auto j = 0U; j < cdepth; j++) {
                        dseq[channels * j + channel] = 0;
                      }
                    }
                }

              T = v_shift_left8(T);
            }

          if (done == sequences) {
            break;
          }

#ifdef __x86_64__
#ifdef __SSE3__
          if (ssse3_present != 0)
            {
              dprofile_shuffle8(dprofile, score_matrix, dseq.data());
            }
          else
#endif
#endif
            {
              dprofile_fill8(dprofile, score_matrix, dseq.data());
            }

          MQ = v_and8(M, Q);
          MR = v_and8(M, R);
          MQ0 = MQ;

          align_cells_masked_8(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
        }

      F0 = v_add8(F0, R);
      F0 = v_add8(F0, R);
      F0 = v_add8(F0, R);
      H0 = v_sub8(F0, Q);
      F0 = v_add8(F0, R);

      dir += 4 * longestdbsequence;
      if (dir >= dirbuffer + dirbuffersize) {
        dir -= dirbuffersize;
      }
    }
}
