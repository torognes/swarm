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

#include "db.h"
#include "utils/backtrack.h"
#include "utils/cpu_features.h"
#include "utils/dseq_fill.h"
#include <array>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t, uint8_t
#include <iterator> // std::next, std::distance
#include <limits>
#include <vector>


// refactoring: C++26 std::simd
#ifdef __aarch64__

#include <arm_neon.h>
#include "arch/aarch64/intrinsics_to_functions.h"
#include "arch/aarch64/search_dispatch.h"
using VECTORTYPE = uint16x8_t;

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
using VECTORTYPE = vector unsigned short;

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
constexpr unsigned int cdepth {4};
constexpr uint8_t n_bits {16};
using BYTE = unsigned char;
using WORD = unsigned short;  // refactoring: uint16_t?

auto dprofile_fill16(WORD * dprofile_word,
                            WORD const * score_matrix,
                            BYTE const * dseq) -> void
{
  static constexpr auto multiplier = 5U;
  static constexpr auto pos0 = 0;
  static constexpr auto pos1 = pos0 + 1;
  static constexpr auto pos2 = pos1 + 1;
  static constexpr auto pos3 = pos2 + 1;
  static constexpr auto pos4 = pos3 + 1;
  static constexpr auto pos5 = pos4 + 1;
  static constexpr auto pos6 = pos5 + 1;
  static constexpr auto pos7 = pos6 + 1;
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
  assert((pos7 * cdepth * channels) + (channels * cdepth) <= max_ptrdiff);
  for (auto j = 0LL; j < cdepth; ++j)
    {
      std::array<unsigned int, channels> d {{}};   // refactoring: name?
      for (auto z = 0U; z < channels; ++z) {
        d[z] = (static_cast<unsigned int>(*std::next(dseq, (j * channels) + z))) << multiplier;
      }

      reg0  = v_load16(cast_vector16(std::next(score_matrix, d[pos0])));
      reg1  = v_load16(cast_vector16(std::next(score_matrix, d[pos1])));
      reg2  = v_load16(cast_vector16(std::next(score_matrix, d[pos2])));
      reg3  = v_load16(cast_vector16(std::next(score_matrix, d[pos3])));
      reg4  = v_load16(cast_vector16(std::next(score_matrix, d[pos4])));
      reg5  = v_load16(cast_vector16(std::next(score_matrix, d[pos5])));
      reg6  = v_load16(cast_vector16(std::next(score_matrix, d[pos6])));
      reg7  = v_load16(cast_vector16(std::next(score_matrix, d[pos7])));

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

      v_store16(cast_vector16(std::next(dprofile_word, (pos0 * cdepth * channels) + (channels * j))), reg24);
      v_store16(cast_vector16(std::next(dprofile_word, (pos1 * cdepth * channels) + (channels * j))), reg25);
      v_store16(cast_vector16(std::next(dprofile_word, (pos2 * cdepth * channels) + (channels * j))), reg26);
      v_store16(cast_vector16(std::next(dprofile_word, (pos3 * cdepth * channels) + (channels * j))), reg27);
      v_store16(cast_vector16(std::next(dprofile_word, (pos4 * cdepth * channels) + (channels * j))), reg28);
      v_store16(cast_vector16(std::next(dprofile_word, (pos5 * cdepth * channels) + (channels * j))), reg29);
      v_store16(cast_vector16(std::next(dprofile_word, (pos6 * cdepth * channels) + (channels * j))), reg30);
      v_store16(cast_vector16(std::next(dprofile_word, (pos7 * cdepth * channels) + (channels * j))), reg31);
    }
}

inline auto onestep_16(VECTORTYPE & H,
                       VECTORTYPE & N,
                       VECTORTYPE & F,
                       VECTORTYPE V,
                       WORD * DIR,
                       VECTORTYPE & E,
                       VECTORTYPE QR,
                       VECTORTYPE R) -> void
{
  H = v_add16(H, V);
  auto W = H;
  H = v_min16(H, F);
  *std::next(DIR, 0) = v_mask_eq16(W, H);
  H = v_min16(H, E);
  *std::next(DIR, 1) = v_mask_eq16(H, E);
  N = H;
  H = v_add16(H, QR);
  F = v_add16(F, R);
  E = v_add16(E, R);
  F = v_min16(H, F);
  *std::next(DIR, 2) = v_mask_eq16(H, F);
  E = v_min16(H, E);
  *std::next(DIR, 3) = v_mask_eq16(H, E);
}


auto align_cells_regular_16(VECTORTYPE * Sm,
                            VECTORTYPE * hep,
                            VECTORTYPE ** qp,
                            VECTORTYPE const * Qm,
                            VECTORTYPE const * Rm,
                            uint64_t ql,
                            VECTORTYPE const * F0,
                            uint64_t * dir_long,
                            VECTORTYPE const * H0) -> void
{
  static constexpr auto step = 16;
  static constexpr auto offset0 = 0;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * dir = reinterpret_cast<WORD *>(dir_long);

  const auto Q = *Qm;
  const auto R = *Rm;

  auto f0 = *F0;
  auto f1 = v_add16(f0, R);
  auto f2 = v_add16(f1, R);
  auto f3 = v_add16(f2, R);

  auto h0 = *H0;
  auto h1 = v_sub16(f0, Q);
  auto h2 = v_add16(h1, R);
  auto h3 = v_add16(h2, R);

  auto h5 = v_zero16();
  auto h6 = v_zero16();
  auto h7 = v_zero16();
  auto h8 = v_zero16();

  assert(ql <= max_ptrdiff);
  assert(ql <= ((max_ptrdiff - 1) / 2));  // max 'E' offset
  assert(ql <= ((max_ptrdiff - offset3) / step));  // max 'dir' offset
  auto const ql_signed = static_cast<std::ptrdiff_t>(ql);
  for (auto pos = 0LL; pos < ql_signed; ++pos)
    {
      VECTORTYPE const * x = *std::next(qp, pos + 0);
      h4 = *std::next(hep, (2 * pos) + 0);
      E  = *std::next(hep, (2 * pos) + 1);
      onestep_16(h0, h5, f0, *std::next(x, 0), std::next(dir, (step * pos) + offset0), E, Q, R);
      onestep_16(h1, h6, f1, *std::next(x, 1), std::next(dir, (step * pos) + offset1), E, Q, R);
      onestep_16(h2, h7, f2, *std::next(x, 2), std::next(dir, (step * pos) + offset2), E, Q, R);
      onestep_16(h3, h8, f3, *std::next(x, 3), std::next(dir, (step * pos) + offset3), E, Q, R);
      *std::next(hep, (2 * pos) + 0) = h8;
      *std::next(hep, (2 * pos) + 1) = E;
      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  *std::next(Sm, 0) = h5;
  *std::next(Sm, 1) = h6;
  *std::next(Sm, 2) = h7;
  *std::next(Sm, 3) = h8;
}


auto align_cells_masked_16(VECTORTYPE * Sm,
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

  auto * dir = reinterpret_cast<WORD *>(dir_long);

  const auto Q = *Qm;
  const auto R = *Rm;

  auto f0 = *F0;
  auto f1 = v_add16(f0, R);
  auto f2 = v_add16(f1, R);
  auto f3 = v_add16(f2, R);

  auto h0 = *H0;
  auto h1 = v_sub16(f0, Q);
  auto h2 = v_add16(h1, R);
  auto h3 = v_add16(h2, R);

  auto h5 = v_zero16();
  auto h6 = v_zero16();
  auto h7 = v_zero16();
  auto h8 = v_zero16();

  assert(ql <= max_ptrdiff);
  assert(ql <= ((max_ptrdiff - 1) / 2));  // max 'E' offset
  assert(ql <= ((max_ptrdiff - offset3) / step));  // max 'dir' offset
  auto const ql_signed = static_cast<std::ptrdiff_t>(ql);
  for (auto pos = 0LL; pos < ql_signed; ++pos)
    {
      h4 = *std::next(hep, (2 * pos) + 0);
      E  = *std::next(hep, (2 * pos) + 1);
      VECTORTYPE const * x = *std::next(qp, pos + 0);

      /* mask h4 and E */
      h4 = v_sub16(h4, *Mm);
      E  = v_sub16(E,  *Mm);

      /* init h4 and E */
      h4 = v_add16(h4, *MQ);
      E  = v_add16(E,  *MQ);
      E  = v_add16(E,  *MQ0);

      /* update MQ */
      *MQ = v_add16(*MQ,  *MR);

      onestep_16(h0, h5, f0, *std::next(x, 0), std::next(dir, (step * pos) + offset0), E, Q, R);
      onestep_16(h1, h6, f1, *std::next(x, 1), std::next(dir, (step * pos) + offset1), E, Q, R);
      onestep_16(h2, h7, f2, *std::next(x, 2), std::next(dir, (step * pos) + offset2), E, Q, R);
      onestep_16(h3, h8, f3, *std::next(x, 3), std::next(dir, (step * pos) + offset3), E, Q, R);
      *std::next(hep, (2 * pos) + 0) = h8;
      *std::next(hep, (2 * pos) + 1) = E;

      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  *std::next(Sm, 0) = h5;
  *std::next(Sm, 1) = h6;
  *std::next(Sm, 2) = h7;
  *std::next(Sm, 3) = h8;
}


namespace {

// Store the final score for the sequence that just ended in 'channel'
// and, when the score fits in a WORD, recover its number of differences
// by backtracking the alignment.
auto save_score_16(int64_t const cand_id,
                          unsigned int const channel,
                          VECTORTYPE const * const score_vectors,
                          std::array<char const *, channels> const & d_address,
                          std::array<uint64_t, channels> const & d_offset,
                          std::array<uint64_t, channels> const & d_length,
                          char const * const qseq,
                          uint64_t const qlen,
                          std::vector<uint64_t> & dirbuffer,
                          uint64_t const q_start_size,
                          uint64_t * const scores,
                          uint64_t * const diffs,
                          uint64_t * const alignmentlengths,
                          uint64_t & done) -> void
{
  static constexpr auto uint16_max = std::numeric_limits<uint16_t>::max();

  // save score

  const uint64_t dbseqlen = d_length[channel];
  const uint64_t z = (dbseqlen + 3) % 4;
  assert(z * channels + channel <= max_ptrdiff);
  const uint64_t score
    = *std::next(reinterpret_cast<WORD const *>(score_vectors), static_cast<std::ptrdiff_t>((z * channels) + channel));
  *std::next(scores, cand_id) = score;

  uint64_t diff {0};

  if (score < uint16_max)
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
      diff = uint16_max;
    }

  *std::next(diffs, cand_id) = diff;

  ++done;
}


// Attach the next database sequence to 'channel': record its address and
// length, reset the per-channel cursors, seed the H0/F0 lanes, and prime
// the first block. Returns whether the channel already reached the end of
// its (short) sequence, i.e. the next block is no longer "easy".
template <std::size_t capacity>
auto load_next_sequence_16(unsigned int const channel,
                                  Data const & data,
                                  uint64_t const * const seqnos,
                                  uint64_t & next_id,
                                  uint64_t const * const dirbuffer_begin,
                                  uint64_t const * const dir,
                                  WORD const gap_open_penalty,
                                  WORD const gap_extend_penalty,
                                  VECTORTYPE & H0,
                                  VECTORTYPE & F0,
                                  std::array<unsigned char, capacity> & dseq,
                                  std::array<int64_t, channels> & seq_id,
                                  std::array<char const *, channels> & d_address,
                                  std::array<uint64_t, channels> & d_length,
                                  std::array<uint64_t, channels> & d_pos,
                                  std::array<uint64_t, channels> & d_offset) -> bool
{
  assert(next_id <= std::numeric_limits<int64_t>::max());
  assert(next_id <= max_ptrdiff);
  // get next sequence
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

  assert(((2U * gap_open_penalty) + (2U * gap_extend_penalty)) <= std::numeric_limits<WORD>::max());
  *std::next(reinterpret_cast<WORD *>(&H0), channel) = 0;
  *std::next(reinterpret_cast<WORD *>(&F0), channel) = static_cast<WORD>((2U * gap_open_penalty) + (2U * gap_extend_penalty));

  // fill channel
  return fill_channel<channels, cdepth>(dseq, channel, d_address, d_pos, d_length);
}

}  // namespace


// search16 is an inherent streaming state machine: in a single pass it
// fills the channels, dispatches the vectorised kernels, and swaps out
// finished database sequences one channel at a time. Extracting the
// fill, save-score and load-next-sequence steps already cut its cognitive
// complexity from 95 to 39; the residual nesting is the per-channel
// switch itself, which is intrinsic to the single-pass design.
auto search16(Data const & data,
              std::vector<WORD *> & q_start,
              WORD gap_open_penalty,
              WORD gap_extend_penalty,
              WORD const * score_matrix,
              std::vector<WORD> & dprofile,
              WORD * hearray,
              uint64_t sequences,
              uint64_t const * seqnos,
              uint64_t * scores,
              uint64_t * diffs,
              uint64_t * alignmentlengths,
              char const * qseq,
              uint64_t qlen,
              std::vector<uint64_t> & dirbuffer,
              Cpu_features const & cpu_features) -> void
{
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
  // an array of BYTES (or WORDS?)
  std::array<BYTE, channels * sizeof(VECTORTYPE) / sizeof(BYTE)> dseq {{}};

  uint64_t next_id {0};
  uint64_t done {0};

  const auto T0 = make_T0_16();

  assert((gap_open_penalty + gap_extend_penalty) <= std::numeric_limits<short>::max());
  assert(gap_extend_penalty <= std::numeric_limits<short>::max());
  auto Q = v_dup16(static_cast<short>(gap_open_penalty + gap_extend_penalty));
  auto R = v_dup16(static_cast<short>(gap_extend_penalty));

  auto * hep = reinterpret_cast<VECTORTYPE *>(hearray);
  auto * * qp = reinterpret_cast<VECTORTYPE * *>(q_start.data());

  auto F0 = v_zero16();
  auto H0 = v_zero16();

  bool easy {false};

  uint64_t * dir = dirbuffer.data();

  while (true) {

      if (easy) {
          // fill all channels

          easy = fill_all_channels<channels, cdepth>(dseq, d_address, d_pos, d_length);

          dispatch_dprofile16(cpu_features, dprofile.data(), score_matrix, dseq.data());

          dispatch_align_regular_16(cpu_features, S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
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
              if (d_pos[channel] < d_length[channel])
                {
                  // this channel has more sequence

                  if (fill_channel<channels, cdepth>(dseq, channel, d_address, d_pos, d_length)) {
                    easy = false;
                  }
                }
              else
                {
                  // sequence in channel ended,
                  // change of sequence

                  M = v_xor16(M, T);

                  const int64_t cand_id = seq_id[channel];

                  if (cand_id >= 0)
                    {
                      save_score_16(cand_id, channel, S,
                                    d_address, d_offset, d_length,
                                    qseq, qlen, dirbuffer, q_start.size(),
                                    scores, diffs, alignmentlengths, done);
                    }

                  if (next_id < sequences)
                    {
                      if (load_next_sequence_16(channel, data, seqnos, next_id,
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

              T = v_shift_left16(T);
            }

          if (done == sequences) {
            break;
          }

          dispatch_dprofile16(cpu_features, dprofile.data(), score_matrix, dseq.data());

          MQ = v_and16(M, Q);
          MR = v_and16(M, R);
          MQ0 = MQ;

          dispatch_align_masked_16(cpu_features, S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
        }

      F0 = v_add16(F0, R);
      F0 = v_add16(F0, R);
      F0 = v_add16(F0, R);
      H0 = v_sub16(F0, Q);
      F0 = v_add16(F0, R);

      assert(4 * q_start.size() <= max_ptrdiff);
      dir = std::next(dir, static_cast<std::ptrdiff_t>(4 * q_start.size()));
      assert(dirbuffer.size() <= max_ptrdiff);
      if (dir >= std::next(dirbuffer.data(), static_cast<std::ptrdiff_t>(dirbuffer.size()))) {
        dir = std::prev(dir, static_cast<std::ptrdiff_t>(dirbuffer.size()));
      }
    }
}
