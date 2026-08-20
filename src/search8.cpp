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

#include "search8.hpp"
#include "db.hpp"
#include "utils/backtrack.hpp"
#include "utils/search_data.hpp"  // Search_data, BYTE, score_ceiling_8 (pulls in Cpu_features)
#include "utils/dseq_fill.hpp"
#include "utils/mask_vectors.hpp"  // No_mask, Mask_vectors
#include "utils/span.hpp"  // Span<uint64_t>
#include "utils/view.hpp"  // View<uint64_t>, make_view
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t, uint8_t
#include <cstring>  // std::memcpy
#include <iterator> // std::next
#include <limits>


// refactoring: C++26 std::simd
#ifdef __aarch64__

#include <arm_neon.h>
#include "arch/aarch64/intrinsics_to_functions.hpp"
#include "arch/aarch64/search_dispatch.hpp"
using VECTORTYPE = uint8x16_t;

#elif defined __x86_64__

#ifdef __SSE2__

#include <emmintrin.h>  // SSE2 intrinsics
#include "arch/x86_64/intrinsics_to_functions.hpp"
using VECTORTYPE = __m128i;

#endif

#include "arch/x86_64/search_dispatch.hpp"

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__

#include <altivec.h>
#include "arch/ppc/intrinsics_to_functions.hpp"
#include "arch/ppc/search_dispatch.hpp"
using VECTORTYPE = __vector unsigned char;

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
static_assert(channels == channels_at_8_bits, "Dseq_8 is sized for 16 channels");
constexpr unsigned int cdepth {4};
static_assert(cdepth == depth_slots, "a block covers exactly one packed byte");
constexpr uint8_t n_bits {bits8};  // backtrack.hpp

// One pass over eight score-matrix columns, where this used to make four
// passes over all thirty-two. The kernel only ever dereferences the profile
// slots for nucleotide codes 1 to 4 -- Scanner::init points qtable_v at
// 64 * code (utils/scanner.cpp) and align_cells reads four vectors from
// there -- so the slots this filled for codes 8 to 31 were written and
// never read. Poisoning the score matrix outside its 4x4 ACGT core leaves
// every output byte-identical; see utils/score_matrix.hpp for what that
// core is. dprofile_fill16 has always written one chunk only, and the two
// builders now have the same shape.
//
// Dropping the other three chunks also made this small enough to inline:
// it was 1,854 bytes of its own before, which is what an earlier note here
// recorded as objdump showing it was not inlined.
//
// The three buffers arrive as their own containers rather than as three
// same-family pointers in a row: their types now say which is which, and
// the '_a'/'_v' parameters are unpacked once here into the base pointers
// the body walks. Inlining .data() at each cast site instead cost two
// prologue instructions and 16 bytes of search8 when the same thing was
// tried on hearray (commit 80df510).
auto dprofile_fill8(Dprofile_8 & dprofile_a,
                           Score_matrix_8 const & score_matrix_a,
                           Dseq_8 const & dseq_a) -> void
{
  auto * const dprofile = dprofile_a.data();
  auto const * const score_matrix = score_matrix_a.data();

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
        score_offsets[i] = (static_cast<unsigned int>(dseq_a[(j * channels) + i])) << multiplier;
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
    }
}


namespace {

// The lane operations and direction-word type for this file's width, handed
// to the shared kernel (utils/align_cells.hpp). Dir_word is 16 bits even at
// this width: v_mask_eq8 packs one bit per channel and there are sixteen.
struct Ops_8 {
  using Dir_word = unsigned short;
  static auto add(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_add8(lhs, rhs); }
  static auto sub(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_sub8(lhs, rhs); }
  static auto min(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_min8(lhs, rhs); }
  static auto mask_eq(VECTORTYPE const lhs, VECTORTYPE const rhs) -> Dir_word { return v_mask_eq8(lhs, rhs); }
  static auto zero() -> VECTORTYPE { return v_zero8(); }
};


// Last, inside this anonymous namespace: the shared kernel names VECTORTYPE,
// No_mask, Ops_8 and max_ptrdiff, all of which are in scope only from here,
// and being inside the namespace gives its instantiations the internal linkage
// the hand-written copy had. See the header for why it cannot include what it
// uses.
#include "utils/align_cells.hpp"

}  // namespace


namespace {

auto align_cells_regular_8(VECTORTYPE * const Sm,
                           VECTORTYPE * const hep,
                           VECTORTYPE ** const qp,
                           VECTORTYPE const & Qm,
                           VECTORTYPE const & Rm,
                           uint64_t const ql,
                           VECTORTYPE const & F0,
                           uint64_t * const dir_long,
                           VECTORTYPE const & H0) -> void
{
  No_mask no_mask;
  align_cells<Ops_8>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, no_mask);
}


auto align_cells_masked_8(VECTORTYPE * const Sm,
                          VECTORTYPE * const hep,
                          VECTORTYPE ** const qp,
                          VECTORTYPE const & Qm,
                          VECTORTYPE const & Rm,
                          uint64_t const ql,
                          VECTORTYPE const & F0,
                          uint64_t * const dir_long,
                          VECTORTYPE const & H0,
                          VECTORTYPE const * const Mm,
                          VECTORTYPE * const MQ,
                          VECTORTYPE const * const MR,
                          VECTORTYPE const * const MQ0) -> void
{
  Mask_vectors masks {*Mm, *MQ, *MR, *MQ0};
  align_cells<Ops_8>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, masks);
  *MQ = masks.mq;
}

}  // namespace


namespace {

// Names the three shared channel operations need, then the operations
// themselves. LANE and score_ceiling are this file's width; VECTORTYPE,
// channels, cdepth, n_bits and max_ptrdiff are already above. The header
// includes nothing and has to come last, inside this namespace -- see it
// for why.
using LANE = BYTE;
constexpr auto score_ceiling = score_ceiling_8;
#include "utils/search_channel_ops.hpp"

}  // namespace


// search8 is an inherent streaming state machine: in a single pass it
// fills the channels, dispatches the vectorised kernels, and swaps out
// finished database sequences one channel at a time. Extracting the
// fill, save-score and load-next-sequence steps already cut its cognitive
// complexity from 95 to 39; the residual nesting is the per-channel
// switch itself, which is intrinsic to the single-pass design.
auto search8(Data const & data,
             Search_data & search_data,
             BYTE const gap_open_penalty,
             BYTE const gap_extend_penalty,
             Score_matrix_8 const & score_matrix,
             View<uint64_t> const seqnos,
             Span<uint64_t> const scores,
             Span<uint64_t> const diffs,
             Sequence const & query) -> void
{
  // the three target arrays are one window over the same candidate list,
  // so they all carry its length; search_data.target_count is that same
  // count, kept for the loop below to read as a number
  assert(scores.size() == seqnos.size());
  assert(diffs.size() == seqnos.size());

  // unpack the per-thread working set (see utils/search_data.hpp)
  auto & q_start = search_data.qtable_v;
  auto & dprofile = search_data.dprofile_a;
  auto * const hearray = search_data.hearray_v.data();  // He_block *
  auto const sequences = seqnos.size();
  auto const qlen = static_cast<uint64_t>(query.length);
  auto const dirbuffer = make_span(search_data.dir_array_v);
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
  std::array<Sequence, channels> d_sequence {{}};
  std::array<int64_t, channels> seq_id {{}};
  seq_id.fill(-1);

  // refactoring fail: std::array -> warning: ignoring attributes on
  // template argument ‘VECTORTYPE’ {aka ‘__m128i’}
  VECTORTYPE S[4];

  // Dseq_8 is an array of size VECTORTYPE * channels, interpreted as
  // an array of BYTES (see utils/search_data.hpp)
  Dseq_8 dseq {};

  uint64_t next_id {0};
  uint64_t done {0};

  auto const T0 = make_T0_8();

  assert(gap_open_penalty + gap_extend_penalty <= std::numeric_limits<char>::max());
  assert(gap_extend_penalty <= std::numeric_limits<char>::max());
  auto Q = v_dup8(static_cast<char>(gap_open_penalty + gap_extend_penalty));
  auto R = v_dup8(static_cast<char>(gap_extend_penalty));

  // refactoring: can't remove reinterpret_cast, cast_vector8() is a nullop in Aarch64
  // hearray is a He_block * now, so this cast no longer widens the alignment
  // of a bare BYTE * -- the source type already carries the 16 the load needs.
  auto *hep = reinterpret_cast<VECTORTYPE*>(hearray);
  auto ** const qp = reinterpret_cast<VECTORTYPE**>(q_start.data());

  auto F0 = v_zero8();
  auto H0 = v_zero8();

  bool easy {false};

  uint64_t * dir = dirbuffer.begin();

  while (true) {
      if (easy) {
          // fill all channels

          easy = fill_all_channels<channels, cdepth>(dseq, d_sequence, d_pos);

          dispatch_dprofile8(cpu_features, dprofile, score_matrix, dseq);

          align_cells_regular_8(S, hep, qp, Q, R, qlen, F0, dir, H0);
        }
      else
        {
          // One or more sequences ended in the previous block
          // We have to switch over to a new sequence

          easy = true;

          M = v_zero8();
          T = T0;
          for (auto channel = 0U; channel < channels; ++channel) {
              if (d_pos[channel] < d_sequence[channel].length) {
                  // this channel has more sequence

                  if (fill_channel<channels, cdepth>(dseq, channel, d_sequence, d_pos)) {
                    easy = false;
                  }
                }
              else
                {
                  // sequence in channel ended,
                  // change of sequence

                  M = v_xor8(M, T);

                  int64_t const cand_id = seq_id[channel];

                  if (cand_id >= 0)
                    {
                      save_score(cand_id, channel, S,
                                   d_sequence, d_offset,
                                   query, static_cast<View<uint64_t>>(dirbuffer), q_start.size(),
                                   scores, diffs, done);
                    }

                  if (next_id < sequences)
                    {
                      if (load_next_sequence(channel, data, seqnos, next_id,
                                               dirbuffer, dir,
                                               gap_open_penalty, gap_extend_penalty,
                                               H0, F0, dseq,
                                               seq_id, d_sequence, d_pos, d_offset)) {
                        easy = false;
                      }
                    }
                  else
                    {
                      // no more sequences, empty channel
                      seq_id[channel] = -1;
                      d_sequence[channel] = Sequence{};
                      d_pos[channel] = 0;
                      clear_channel<channels, cdepth>(dseq, channel);
                    }
                }

              T = v_shift_left8(T);
            }

          if (done == sequences) {
            break;
          }

          dispatch_dprofile8(cpu_features, dprofile, score_matrix, dseq);

          MQ = v_and8(M, Q);
          MR = v_and8(M, R);
          MQ0 = MQ;

          align_cells_masked_8(S, hep, qp, Q, R, qlen, F0, dir, H0, &M, &MQ, &MR, &MQ0);
        }

      F0 = v_add8(F0, R);
      F0 = v_add8(F0, R);
      F0 = v_add8(F0, R);
      H0 = v_sub8(F0, Q);
      F0 = v_add8(F0, R);

      assert(4 * q_start.size() <= max_ptrdiff);
      dir = std::next(dir, static_cast<std::ptrdiff_t>(4 * q_start.size()));
      // the cursor sweeps the buffer and wraps at its end; both bounds
      // come from the span rather than from pointer arithmetic on .data()
      if (dir >= dirbuffer.end()) {
        dir = std::prev(dir, static_cast<std::ptrdiff_t>(dirbuffer.size()));
      }
    }
}
