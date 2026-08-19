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

#include "search16.hpp"
#include "db.hpp"
#include "utils/backtrack.hpp"
#include "utils/search_data.hpp"  // Search_data, WORD, score_ceiling_16 (pulls in Cpu_features)
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
using VECTORTYPE = uint16x8_t;

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
using VECTORTYPE = __vector unsigned short;

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

constexpr unsigned int channels {8};
static_assert(channels == channels_at_16_bits, "Dseq_16 is sized for 8 channels");
constexpr unsigned int cdepth {4};
static_assert(cdepth == depth_slots, "a block covers exactly one packed byte");
constexpr uint8_t n_bits {bits16};  // backtrack.hpp

// The three buffers arrive as their own containers rather than as three
// same-family pointers in a row: their types now say which is which, and
// the '_a'/'_v' parameters are unpacked once here into the base pointers
// the body walks. Inlining .data() at each cast site instead cost two
// prologue instructions and 16 bytes of search8 when the same thing was
// tried on hearray (commit 80df510).
auto dprofile_fill16(Dprofile_16 & dprofile_a,
                            Score_matrix_16 const & score_matrix_a,
                            Dseq_16 const & dseq_a) -> void
{
  auto * const dprofile_word = dprofile_a.data();
  auto const * const score_matrix = score_matrix_a.data();

  static constexpr auto multiplier = 5U;
  static constexpr auto pos0 = 0;
  static constexpr auto pos1 = pos0 + 1;
  static constexpr auto pos2 = pos1 + 1;
  static constexpr auto pos3 = pos2 + 1;
  static constexpr auto pos4 = pos3 + 1;
  static constexpr auto pos5 = pos4 + 1;
  static constexpr auto pos6 = pos5 + 1;
  static constexpr auto pos7 = pos6 + 1;
  static constexpr auto offset0 = pos0 * cdepth * channels;
  static constexpr auto offset1 = pos1 * cdepth * channels;
  static constexpr auto offset2 = pos2 * cdepth * channels;
  static constexpr auto offset3 = pos3 * cdepth * channels;
  static constexpr auto offset4 = pos4 * cdepth * channels;
  static constexpr auto offset5 = pos5 * cdepth * channels;
  static constexpr auto offset6 = pos6 * cdepth * channels;
  static constexpr auto offset7 = pos7 * cdepth * channels;
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

  assert(cdepth <= ((max_ptrdiff - channels) / channels));  // max 'd' offset
  assert(channels <= std::numeric_limits<long int>::max());
  assert(channels <= std::numeric_limits<unsigned int>::max());
  assert(offset7 + (static_cast<std::ptrdiff_t>(channels) * cdepth) <= max_ptrdiff);
  for (auto j = 0LL; j < cdepth; ++j)
    {
      // dseq is laid out as cdepth blocks of 'channels' bytes; j is a
      // long long here because 'lane' below is a std::ptrdiff_t
      auto const block = static_cast<std::size_t>(j) * channels;
      std::array<unsigned int, channels> score_offsets {{}};
      for (auto z = 0U; z < channels; ++z) {
        score_offsets[z] = (static_cast<unsigned int>(dseq_a[block + z])) << multiplier;
      }

      reg0  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos0])));
      reg1  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos1])));
      reg2  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos2])));
      reg3  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos3])));
      reg4  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos4])));
      reg5  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos5])));
      reg6  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos6])));
      reg7  = v_load16(cast_vector16(std::next(score_matrix, score_offsets[pos7])));

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

      std::ptrdiff_t const lane = static_cast<std::ptrdiff_t>(channels) * j;
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset0)), reg24);
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset1)), reg25);
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset2)), reg26);
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset3)), reg27);
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset4)), reg28);
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset5)), reg29);
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset6)), reg30);
      v_store16(cast_vector16(std::next(dprofile_word, lane + offset7)), reg31);
    }
}

namespace {

// The lane operations and direction-word type for this file's width, handed
// to the shared kernel (utils/align_cells.hpp). min() is whatever v_min16 is on this architecture; on x86-64
// without SSE4.1 that is an emulation, which is why arch/x86_64/sse41.cpp
// instantiates the same kernel with PMINUW instead.
struct Ops_baseline {
  using Dir_word = WORD;
  static auto add(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_add16(lhs, rhs); }
  static auto sub(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_sub16(lhs, rhs); }
  static auto min(VECTORTYPE const lhs, VECTORTYPE const rhs) -> VECTORTYPE { return v_min16(lhs, rhs); }
  static auto mask_eq(VECTORTYPE const lhs, VECTORTYPE const rhs) -> Dir_word { return v_mask_eq16(lhs, rhs); }
  static auto zero() -> VECTORTYPE { return v_zero16(); }
};




// Last, inside this anonymous namespace: the shared kernel names
// VECTORTYPE, WORD, the v_* wrappers, Mask_vectors, apply_mask and
// max_ptrdiff, all of which are in scope only from here, and being inside
// the namespace gives its instantiations the internal linkage the two
// hand-written copies had. See the header for why it cannot include what
// it uses.
#include "utils/align_cells.hpp"

}  // namespace


auto align_cells_regular_16(VECTORTYPE * const Sm,
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
  align_cells<Ops_baseline>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, no_mask);
}


auto align_cells_masked_16(VECTORTYPE * const Sm,
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
  align_cells<Ops_baseline>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, masks);
  *MQ = masks.mq;
}


namespace {

// Names the three shared channel operations need, then the operations
// themselves. LANE and score_ceiling are this file's width; VECTORTYPE,
// channels, cdepth, n_bits and max_ptrdiff are already above. The header
// includes nothing and has to come last, inside this namespace -- see it
// for why.
using LANE = WORD;
constexpr auto score_ceiling = score_ceiling_16;
#include "utils/search_channel_ops.hpp"

}  // namespace


// search16 is an inherent streaming state machine: in a single pass it
// fills the channels, dispatches the vectorised kernels, and swaps out
// finished database sequences one channel at a time. Extracting the
// fill, save-score and load-next-sequence steps already cut its cognitive
// complexity from 95 to 39; the residual nesting is the per-channel
// switch itself, which is intrinsic to the single-pass design.
auto search16(Data const & data,
              Search_data & search_data,
              WORD const gap_open_penalty,
              WORD const gap_extend_penalty,
              Score_matrix_16 const & score_matrix,
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
  auto & q_start = search_data.qtable_w_v;
  auto & dprofile = search_data.dprofile_w_a;
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

  // Dseq_16 is an array of size VECTORTYPE * channels, interpreted as
  // an array of BYTES (or WORDS?) -- see utils/search_data.hpp
  Dseq_16 dseq {};

  uint64_t next_id {0};
  uint64_t done {0};

  auto const T0 = make_T0_16();

  assert((gap_open_penalty + gap_extend_penalty) <= std::numeric_limits<short>::max());
  assert(gap_extend_penalty <= std::numeric_limits<short>::max());
  auto Q = v_dup16(static_cast<short>(gap_open_penalty + gap_extend_penalty));
  auto R = v_dup16(static_cast<short>(gap_extend_penalty));

  // one cast, from the over-aligned He_block straight to the vector type.
  // hearray used to be a WORD * that nothing else read, and that step
  // discarded the alignment these loads need.
  auto * hep = reinterpret_cast<VECTORTYPE *>(hearray);
  auto * * const qp = reinterpret_cast<VECTORTYPE * *>(q_start.data());

  auto F0 = v_zero16();
  auto H0 = v_zero16();

  bool easy {false};

  uint64_t * dir = dirbuffer.begin();

  while (true) {

      if (easy) {
          // fill all channels

          easy = fill_all_channels<channels, cdepth>(dseq, d_sequence, d_pos);

          dispatch_dprofile16(cpu_features, dprofile, score_matrix, dseq);

          dispatch_align_regular_16(cpu_features, S, hep, qp, Q, R, qlen, F0, dir, H0);
        }
      else
        {
          // One or more sequences ended in the previous block
          // We have to switch over to a new sequence

          easy = true;

          M = v_zero16();
          T = T0;
          for (auto channel = 0U; channel < channels; ++channel)
            {
              if (d_pos[channel] < d_sequence[channel].length)
                {
                  // this channel has more sequence

                  if (fill_channel<channels, cdepth>(dseq, channel, d_sequence, d_pos)) {
                    easy = false;
                  }
                }
              else
                {
                  // sequence in channel ended,
                  // change of sequence

                  M = v_xor16(M, T);

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

              T = v_shift_left16(T);
            }

          if (done == sequences) {
            break;
          }

          dispatch_dprofile16(cpu_features, dprofile, score_matrix, dseq);

          MQ = v_and16(M, Q);
          MR = v_and16(M, R);
          MQ0 = MQ;

          dispatch_align_masked_16(cpu_features, S, hep, qp, Q, R, qlen, F0, dir, H0, &M, &MQ, &MR, &MQ0);
        }

      F0 = v_add16(F0, R);
      F0 = v_add16(F0, R);
      F0 = v_add16(F0, R);
      H0 = v_sub16(F0, Q);
      F0 = v_add16(F0, R);

      assert(4 * q_start.size() <= max_ptrdiff);
      dir = std::next(dir, static_cast<std::ptrdiff_t>(4 * q_start.size()));
      // the cursor sweeps the buffer and wraps at its end; both bounds
      // come from the span rather than from pointer arithmetic on .data()
      if (dir >= dirbuffer.end()) {
        dir = std::prev(dir, static_cast<std::ptrdiff_t>(dirbuffer.size()));
      }
    }
}
