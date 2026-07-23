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

#include "search8.h"
#include "db.h"
#include "utils/backtrack.h"
#include "utils/search_data.h"  // Search_data (pulls in Cpu_features)
#include "utils/dseq_fill.h"
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t, uint8_t
#include <iterator> // std::next
#include <limits>
#include <vector>


// refactoring: C++26 std::simd
#ifdef __aarch64__

#include <arm_neon.h>
#include "arch/aarch64/intrinsics_to_functions.h"
#include "arch/aarch64/search_dispatch.h"
using VECTORTYPE = uint8x16_t;

#elif defined __x86_64__

#ifdef __SSE2__

#include <emmintrin.h>  // SSE2 intrinsics
#include "arch/x86_64/intrinsics_to_functions.h"
using VECTORTYPE = __m128i;

#endif

#include "arch/x86_64/search_dispatch.h"

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__

#include <altivec.h>
#include "arch/ppc/intrinsics_to_functions.h"
#include "arch/ppc/search_dispatch.h"
using VECTORTYPE = vector unsigned char;

#else

#error Big endian ppc64 CPUs not supported

#endif

#else

#error Unknown architecture

#endif


#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
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
auto dprofile_fill8(BYTE * dprofile,
                           BYTE const * score_matrix,
                           BYTE const * dseq) -> void
{
  static constexpr auto multiplier = 5U;
  static_assert((std::numeric_limits<BYTE>::max() << multiplier) <= std::numeric_limits<unsigned int>::max(),
                "score-matrix byte offset must fit in an unsigned int");

  static constexpr auto pos0  = 0U;
  static constexpr auto pos1  = pos0  + 1;
  static constexpr auto pos2  = pos1  + 1;
  static constexpr auto pos3  = pos2  + 1;
  static constexpr auto pos4  = pos3  + 1;
  static constexpr auto pos5  = pos4  + 1;
  static constexpr auto pos6  = pos5  + 1;
  static constexpr auto pos7  = pos6  + 1;
  static constexpr auto pos8  = pos7  + 1;
  static constexpr auto pos9  = pos8  + 1;
  static constexpr auto pos10 = pos9  + 1;
  static constexpr auto pos11 = pos10 + 1;
  static constexpr auto pos12 = pos11 + 1;
  static constexpr auto pos13 = pos12 + 1;
  static constexpr auto pos14 = pos13 + 1;
  static constexpr auto pos15 = pos14 + 1;

  static constexpr auto line0  = 64U * 0;  // as in 'cache line': 64 bytes
  static constexpr auto line1  = 64U * 1;
  static constexpr auto line2  = 64U * 2;
  static constexpr auto line3  = 64U * 3;
  static constexpr auto line4  = 64U * 4;
  static constexpr auto line5  = 64U * 5;
  static constexpr auto line6  = 64U * 6;
  static constexpr auto line7  = 64U * 7;
  static constexpr auto line8  = 64U * 8;
  static constexpr auto line16 = 64U * 16;  // 1,024
  static constexpr auto line24 = 64U * 24;  // 1,536

  static constexpr auto offset8  =  8U;
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

  for (auto j = 0U; j < cdepth; ++j)
    {
      std::array<unsigned int, channels> score_offsets {{}};
      for (auto i = 0U; i < channels; ++i) {
        score_offsets[i] = (static_cast<unsigned int>(*std::next(dseq, (j * channels) + i))) << multiplier;
      }

      reg0  = v_load_64(std::next(score_matrix, score_offsets[pos0]));
      reg2  = v_load_64(std::next(score_matrix, score_offsets[pos2]));
      reg4  = v_load_64(std::next(score_matrix, score_offsets[pos4]));
      reg6  = v_load_64(std::next(score_matrix, score_offsets[pos6]));
      reg8  = v_load_64(std::next(score_matrix, score_offsets[pos8]));
      reg10 = v_load_64(std::next(score_matrix, score_offsets[pos10]));
      reg12 = v_load_64(std::next(score_matrix, score_offsets[pos12]));
      reg14 = v_load_64(std::next(score_matrix, score_offsets[pos14]));

      reg0  = v_merge_lo_8(reg0,  *cast_vector8(std::next(score_matrix, score_offsets[pos1])));
      reg2  = v_merge_lo_8(reg2,  *cast_vector8(std::next(score_matrix, score_offsets[pos3])));
      reg4  = v_merge_lo_8(reg4,  *cast_vector8(std::next(score_matrix, score_offsets[pos5])));
      reg6  = v_merge_lo_8(reg6,  *cast_vector8(std::next(score_matrix, score_offsets[pos7])));
      reg8  = v_merge_lo_8(reg8,  *cast_vector8(std::next(score_matrix, score_offsets[pos9])));
      reg10 = v_merge_lo_8(reg10, *cast_vector8(std::next(score_matrix, score_offsets[pos11])));
      reg12 = v_merge_lo_8(reg12, *cast_vector8(std::next(score_matrix, score_offsets[pos13])));
      reg14 = v_merge_lo_8(reg14, *cast_vector8(std::next(score_matrix, score_offsets[pos15])));

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

      std::ptrdiff_t const lane = static_cast<std::ptrdiff_t>(channels) * j;
      v_store8(cast_vector8(std::next(dprofile, lane + line0)), reg0);
      v_store8(cast_vector8(std::next(dprofile, lane + line1)), reg3);
      v_store8(cast_vector8(std::next(dprofile, lane + line2)), reg2);
      v_store8(cast_vector8(std::next(dprofile, lane + line3)), reg7);
      v_store8(cast_vector8(std::next(dprofile, lane + line4)), reg1);
      v_store8(cast_vector8(std::next(dprofile, lane + line5)), reg11);
      v_store8(cast_vector8(std::next(dprofile, lane + line6)), reg6);
      v_store8(cast_vector8(std::next(dprofile, lane + line7)), reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos0]));
      reg1  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos1]));
      reg2  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos2]));
      reg3  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos3]));
      reg4  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos4]));
      reg5  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos5]));
      reg6  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos6]));
      reg7  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos7]));
      reg8  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos8]));
      reg9  = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos9]));
      reg10 = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos10]));
      reg11 = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos11]));
      reg12 = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos12]));
      reg13 = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos13]));
      reg14 = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos14]));
      reg15 = v_load_64(std::next(score_matrix, offset8 + score_offsets[pos15]));

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

      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line0)), reg0);
      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line1)), reg3);
      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line2)), reg2);
      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line3)), reg7);
      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line4)), reg1);
      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line5)), reg11);
      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line6)), reg6);
      v_store8(cast_vector8(std::next(dprofile, lane + line8 + line7)), reg15);


      reg0  = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos0]));
      reg2  = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos2]));
      reg4  = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos4]));
      reg6  = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos6]));
      reg8  = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos8]));
      reg10 = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos10]));
      reg12 = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos12]));
      reg14 = v_load_64(std::next(score_matrix, offset16 + score_offsets[pos14]));

      reg0  = v_merge_lo_8(reg0,  *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos1])));
      reg2  = v_merge_lo_8(reg2,  *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos3])));
      reg4  = v_merge_lo_8(reg4,  *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos5])));
      reg6  = v_merge_lo_8(reg6,  *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos7])));
      reg8  = v_merge_lo_8(reg8,  *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos9])));
      reg10 = v_merge_lo_8(reg10, *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos11])));
      reg12 = v_merge_lo_8(reg12, *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos13])));
      reg14 = v_merge_lo_8(reg14, *cast_vector8(std::next(score_matrix, offset16 + score_offsets[pos15])));

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

      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line0)), reg0);
      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line1)), reg3);
      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line2)), reg2);
      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line3)), reg7);
      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line4)), reg1);
      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line5)), reg11);
      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line6)), reg6);
      v_store8(cast_vector8(std::next(dprofile, lane + line16 + line7)), reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos0]));
      reg1  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos1]));
      reg2  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos2]));
      reg3  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos3]));
      reg4  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos4]));
      reg5  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos5]));
      reg6  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos6]));
      reg7  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos7]));
      reg8  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos8]));
      reg9  = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos9]));
      reg10 = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos10]));
      reg11 = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos11]));
      reg12 = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos12]));
      reg13 = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos13]));
      reg14 = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos14]));
      reg15 = v_load_64(std::next(score_matrix, offset24 + score_offsets[pos15]));

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

      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line0)), reg0);
      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line1)), reg3);
      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line2)), reg2);
      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line3)), reg7);
      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line4)), reg1);
      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line5)), reg11);
      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line6)), reg6);
      v_store8(cast_vector8(std::next(dprofile, lane + line24 + line7)), reg15);
    }
}


namespace {

inline auto onestep_8(VECTORTYPE & H,
                      VECTORTYPE & N,
                      VECTORTYPE & F,
                      VECTORTYPE V,
                      unsigned short * DIR,
                      VECTORTYPE & E,
                      VECTORTYPE QR,
                      VECTORTYPE R) -> void
{
  H = v_add8(H, V);
  const auto W = H;
  H = v_min8(H, F);
  DIR[0] = v_mask_eq8(W, H);  // subscript, not std::next: hot loop, see align_cells
  H = v_min8(H, E);
  DIR[1] = v_mask_eq8(H, E);
  N = H;
  H = v_add8(H, QR);
  F = v_add8(F, R);
  E = v_add8(E, R);
  F = v_min8(H, F);
  DIR[2] = v_mask_eq8(H, F);
  E = v_min8(H, E);
  DIR[3] = v_mask_eq8(H, E);
}


// One block of cells, shared by the regular and masked kernels. The
// masked variant differs only by a per-iteration adjustment of h4 and E
// (consuming the Mm / MQ / MR / MQ0 vectors); 'masked' is a compile-time
// flag, so the regular instantiation drops that block entirely and never
// dereferences the (null) masking pointers.
template <bool masked>
auto align_cells_8(VECTORTYPE * Sm,
                   VECTORTYPE * hep,
                   VECTORTYPE ** qp,
                   VECTORTYPE const * Qm,
                   VECTORTYPE const * Rm,
                   uint64_t ql,
                   VECTORTYPE const * F0,
                   uint64_t * dir_long,
                   VECTORTYPE const * H0,
                   VECTORTYPE const * Mm,
                   VECTORTYPE * MQ,
                   VECTORTYPE const * MR,
                   VECTORTYPE const * MQ0) -> void
{
  static constexpr auto step = 16;
  static constexpr auto offset0 = 0;
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

  assert(ql <= max_ptrdiff);
  assert(ql <= ((max_ptrdiff - 1) / 2));  // max 'E' offset
  assert(ql <= ((max_ptrdiff - offset3) / step));  // max 'dir' offset
  auto const ql_signed = static_cast<std::ptrdiff_t>(ql);
  // Performance: subscript / &dir[...] rather than std::next() in this hot loop
  // (same regression as commit 8c6925f on the 16-bit kernel). Stays clang-tidy
  // clean: pos is signed and operator[] is not pointer arithmetic.
  for (auto pos = 0LL; pos < ql_signed; ++pos)
    {
      VECTORTYPE const * x = qp[pos];
      h4 = hep[(2 * pos) + 0];
      E  = hep[(2 * pos) + 1];

      if (masked)
        {
          /* mask h4 and E */
          h4 = v_sub8(h4, *Mm);
          E  = v_sub8(E,  *Mm);

          /* init h4 and E */
          h4 = v_add8(h4, *MQ);
          E  = v_add8(E,  *MQ);
          E  = v_add8(E,  *MQ0);

          /* update MQ */
          *MQ = v_add8(*MQ,  *MR);
        }

      onestep_8(h0, h5, f0, x[0], &dir[(step * pos) + offset0], E, Q, R);
      onestep_8(h1, h6, f1, x[1], &dir[(step * pos) + offset1], E, Q, R);
      onestep_8(h2, h7, f2, x[2], &dir[(step * pos) + offset2], E, Q, R);
      onestep_8(h3, h8, f3, x[3], &dir[(step * pos) + offset3], E, Q, R);
      hep[(2 * pos) + 0] = h8;
      hep[(2 * pos) + 1] = E;
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

}  // namespace


namespace {

auto align_cells_regular_8(VECTORTYPE * Sm,
                           VECTORTYPE * hep,
                           VECTORTYPE ** qp,
                           VECTORTYPE const * Qm,
                           VECTORTYPE const * Rm,
                           uint64_t ql,
                           VECTORTYPE const * F0,
                           uint64_t * dir_long,
                           VECTORTYPE const * H0) -> void
{
  align_cells_8<false>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0,
                       nullptr, nullptr, nullptr, nullptr);
}


auto align_cells_masked_8(VECTORTYPE * Sm,
                          VECTORTYPE * hep,
                          VECTORTYPE ** qp,
                          VECTORTYPE const * Qm,
                          VECTORTYPE const * Rm,
                          uint64_t ql,
                          VECTORTYPE const * F0,
                          uint64_t * dir_long,
                          VECTORTYPE const * H0,
                          VECTORTYPE const * Mm,
                          VECTORTYPE * MQ,
                          VECTORTYPE const * MR,
                          VECTORTYPE const * MQ0) -> void
{
  align_cells_8<true>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, Mm, MQ, MR, MQ0);
}

}  // namespace


namespace {

// Store the final score for the sequence that just ended in 'channel'
// and, when the score fits in a BYTE, recover its number of differences
// by backtracking the alignment.
auto save_score_8(int64_t const cand_id,
                         unsigned int const channel,
                         VECTORTYPE const * const score_vectors,
                         std::array<char const *, channels> const & d_address,
                         std::array<uint64_t, channels> const & d_offset,
                         std::array<uint64_t, channels> const & d_length,
                         char const * const qseq,
                         uint64_t const qlen,
                         std::vector<uint64_t> const & dirbuffer,
                         uint64_t const q_start_size,
                         uint64_t * const scores,
                         uint64_t * const diffs,
                         uint64_t * const alignmentlengths,
                         uint64_t & done) -> void
{
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();

  // save score

  const uint64_t dbseqlen = d_length[channel];
  const uint64_t z = (dbseqlen + 3) % 4;
  assert(z * channels + channel <= max_ptrdiff);
  const uint64_t score
    = *std::next(reinterpret_cast<BYTE const *>(score_vectors), static_cast<std::ptrdiff_t>((z * channels) + channel));
  *std::next(scores, cand_id) = score;

  uint64_t diff {0};

  if (score < uint8_max)
    {
      char const * dbseq = d_address[channel];
      const uint64_t offset = d_offset[channel];
      diff = backtrack<n_bits>(qseq, dbseq, qlen, dbseqlen,
                               dirbuffer,
                               offset,
                               channel,
                               std::next(alignmentlengths, cand_id),
                               q_start_size);
    }
  else
    {
      diff = uint8_max;
    }

  *std::next(diffs, cand_id) = diff;

  ++done;
}


// Attach the next database sequence to 'channel': record its address and
// length, reset the per-channel cursors, seed the H0/F0 lanes, and prime
// the first block. Returns whether the channel already reached the end of
// its (short) sequence, i.e. the next block is no longer "easy".
template <std::size_t capacity>
auto load_next_sequence_8(unsigned int const channel,
                                 Data const & data,
                                 uint64_t const * const seqnos,
                                 uint64_t & next_id,
                                 uint64_t const * const dirbuffer_begin,
                                 uint64_t const * const dir,
                                 BYTE const gap_open_penalty,
                                 BYTE const gap_extend_penalty,
                                 VECTORTYPE & H0,
                                 VECTORTYPE & F0,
                                 std::array<unsigned char, capacity> & dseq,
                                 std::array<int64_t, channels> & seq_id,
                                 std::array<char const *, channels> & d_address,
                                 std::array<uint64_t, channels> & d_length,
                                 std::array<uint64_t, channels> & d_pos,
                                 std::array<uint64_t, channels> & d_offset) -> bool
{
  // get next sequence
  assert(next_id <= std::numeric_limits<int64_t>::max());
  assert(next_id <= max_ptrdiff);
  seq_id[channel] = static_cast<int64_t>(next_id);
  const uint64_t seqno = *std::next(seqnos, static_cast<std::ptrdiff_t>(next_id));
  auto const & info = data.info(seqno);
  char const * address = info.seq;
  unsigned int const length = info.seqlen;

  d_address[channel] = address;
  d_length[channel] = length;

  d_pos[channel] = 0;
  d_offset[channel] = static_cast<uint64_t>(dir - dirbuffer_begin);
  ++next_id;

  *std::next(reinterpret_cast<BYTE *>(&H0), channel) = 0;
  assert((2U * gap_open_penalty) + (2U * gap_extend_penalty) <= std::numeric_limits<BYTE>::max());
  *std::next(reinterpret_cast<BYTE *>(&F0), channel) = static_cast<BYTE>((2U * gap_open_penalty) + (2U * gap_extend_penalty));

  // fill channel
  return fill_channel<channels, cdepth>(dseq, channel, d_address, d_pos, d_length);
}

}  // namespace


// search8 is an inherent streaming state machine: in a single pass it
// fills the channels, dispatches the vectorised kernels, and swaps out
// finished database sequences one channel at a time. Extracting the
// fill, save-score and load-next-sequence steps already cut its cognitive
// complexity from 95 to 39; the residual nesting is the per-channel
// switch itself, which is intrinsic to the single-pass design.
auto search8(Data const & data,
             Search_data & search_data,
             BYTE gap_open_penalty,
             BYTE gap_extend_penalty,
             BYTE const * score_matrix,
             uint64_t const * seqnos,
             uint64_t * scores,
             uint64_t * diffs,
             uint64_t * alignmentlengths,
             char const * qseq,
             uint64_t qlen) -> void
{
  // unpack the per-thread working set (see utils/search_data.h)
  auto & q_start = search_data.qtable_v;
  auto & dprofile = search_data.dprofile_v;
  auto * const hearray = search_data.hearray_v.data();
  auto const sequences = search_data.target_count;
  auto & dirbuffer = search_data.dir_array_v;
  auto const & cpu_features = search_data.cpu_features;

  VECTORTYPE T;
  VECTORTYPE M;
  VECTORTYPE MQ;
  VECTORTYPE MR;
  VECTORTYPE MQ0;

  // by default, std::array is value-initialized (set to 0 for int,
  // nullptr for pointers, etc)
  std::array<uint64_t, channels> d_pos {{}};
  std::array<uint64_t, channels> d_offset {{}};
  std::array<char const *, channels> d_address {{}};
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

  const auto T0 = make_T0_8();

  assert(gap_open_penalty + gap_extend_penalty <= std::numeric_limits<char>::max());
  assert(gap_extend_penalty <= std::numeric_limits<char>::max());
  auto Q = v_dup8(static_cast<char>(gap_open_penalty + gap_extend_penalty));
  auto R = v_dup8(static_cast<char>(gap_extend_penalty));

  // refactoring: can't remove reinterpret_cast, cast_vector8() is a nullop in Aarch64
  auto *hep = reinterpret_cast<VECTORTYPE*>(hearray);
  auto **qp = reinterpret_cast<VECTORTYPE**>(q_start.data());

  auto F0 = v_zero8();
  auto H0 = v_zero8();

  bool easy {false};

  uint64_t * dir = dirbuffer.data();

  while (true) {
      if (easy) {
          // fill all channels

          easy = fill_all_channels<channels, cdepth>(dseq, d_address, d_pos, d_length);

          dispatch_dprofile8(cpu_features, dprofile.data(), score_matrix, dseq.data());

          align_cells_regular_8(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
        }
      else
        {
          // One or more sequences ended in the previous block
          // We have to switch over to a new sequence

          easy = true;

          M = v_zero8();
          T = T0;
          for (auto channel = 0U; channel < channels; ++channel) {
              if (d_pos[channel] < d_length[channel]) {
                  // this channel has more sequence

                  if (fill_channel<channels, cdepth>(dseq, channel, d_address, d_pos, d_length)) {
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
                      save_score_8(cand_id, channel, S,
                                   d_address, d_offset, d_length,
                                   qseq, qlen, dirbuffer, q_start.size(),
                                   scores, diffs, alignmentlengths, done);
                    }

                  if (next_id < sequences)
                    {
                      if (load_next_sequence_8(channel, data, seqnos, next_id,
                                               dirbuffer.data(), dir,
                                               gap_open_penalty, gap_extend_penalty,
                                               H0, F0, dseq,
                                               seq_id, d_address, d_length, d_pos, d_offset)) {
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
                      clear_channel<channels, cdepth>(dseq, channel);
                    }
                }

              T = v_shift_left8(T);
            }

          if (done == sequences) {
            break;
          }

          dispatch_dprofile8(cpu_features, dprofile.data(), score_matrix, dseq.data());

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

      assert(4 * q_start.size() <= max_ptrdiff);
      dir = std::next(dir, static_cast<std::ptrdiff_t>(4 * q_start.size()));
      assert(dirbuffer.size() <= max_ptrdiff);
      if (dir >= std::next(dirbuffer.data(), static_cast<std::ptrdiff_t>(dirbuffer.size()))) {
        dir = std::prev(dir, static_cast<std::ptrdiff_t>(dirbuffer.size()));
      }
    }
}
