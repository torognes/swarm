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

#include "db.h"
#include "utils/backtrack.h"
#include "utils/queryinfo.h"
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
#include "utils/intrinsics_to_functions_aarch64.h"
using VECTORTYPE = uint16x8_t;

#elif defined __x86_64__

#ifdef __SSE2__

#include <emmintrin.h>  // SSE2 intrinsics
#include "utils/intrinsics_to_functions_x86_64.h"
#include "utils/x86_cpu_feature_ssse3.h"
#include "utils/x86_cpu_feature_sse41.h"
using VECTORTYPE = __m128i;

#endif

#ifdef __SSE3__

#include "ssse3.h"

#endif

#ifdef __SSE4_1__

#include "sse41.h"

#endif

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__

#include <altivec.h>
#include "utils/intrinsics_to_functions_ppc.h"
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

inline auto dprofile_fill16(WORD * dprofile_word,
                            WORD * score_matrix,
                            BYTE const * dseq) -> void
{
  static constexpr auto s_channels = static_cast<int>(channels);
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
  assert((channels + pos7) * cdepth * channels + channels * cdepth <= max_ptrdiff);
  for(auto j = 0LL; j < cdepth; ++j)
    {
      std::array<unsigned int, channels> d {{}};   // refactoring: name?
      for(auto z = 0U; z < channels; ++z) {
        d[z] = (static_cast<unsigned int>(*std::next(dseq, j * channels + z))) << multiplier;
      }

      for(auto i = 0L; i < s_channels; i += s_channels)
        {
          reg0  = v_load16(cast_vector16(std::next(score_matrix, d[pos0] + i)));
          reg1  = v_load16(cast_vector16(std::next(score_matrix, d[pos1] + i)));
          reg2  = v_load16(cast_vector16(std::next(score_matrix, d[pos2] + i)));
          reg3  = v_load16(cast_vector16(std::next(score_matrix, d[pos3] + i)));
          reg4  = v_load16(cast_vector16(std::next(score_matrix, d[pos4] + i)));
          reg5  = v_load16(cast_vector16(std::next(score_matrix, d[pos5] + i)));
          reg6  = v_load16(cast_vector16(std::next(score_matrix, d[pos6] + i)));
          reg7  = v_load16(cast_vector16(std::next(score_matrix, d[pos7] + i)));

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

          v_store16(cast_vector16(std::next(dprofile_word, (i + pos0) * cdepth * channels + channels * j)), reg24);
          v_store16(cast_vector16(std::next(dprofile_word, (i + pos1) * cdepth * channels + channels * j)), reg25);
          v_store16(cast_vector16(std::next(dprofile_word, (i + pos2) * cdepth * channels + channels * j)), reg26);
          v_store16(cast_vector16(std::next(dprofile_word, (i + pos3) * cdepth * channels + channels * j)), reg27);
          v_store16(cast_vector16(std::next(dprofile_word, (i + pos4) * cdepth * channels + channels * j)), reg28);
          v_store16(cast_vector16(std::next(dprofile_word, (i + pos5) * cdepth * channels + channels * j)), reg29);
          v_store16(cast_vector16(std::next(dprofile_word, (i + pos6) * cdepth * channels + channels * j)), reg30);
          v_store16(cast_vector16(std::next(dprofile_word, (i + pos7) * cdepth * channels + channels * j)), reg31);
        }
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
                            VECTORTYPE * Qm,
                            VECTORTYPE * Rm,
                            uint64_t ql,
                            VECTORTYPE * F0,
                            uint64_t * dir_long,
                            VECTORTYPE * H0) -> void
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
  for(auto pos = 0LL; pos < ql_signed; ++pos)
    {
      VECTORTYPE * x = *std::next(qp, pos + 0);
      h4 = *std::next(hep, 2 * pos + 0);
      E  = *std::next(hep, 2 * pos + 1);
      onestep_16(h0, h5, f0, *std::next(x, 0), std::next(dir, step * pos + offset0), E, Q, R);
      onestep_16(h1, h6, f1, *std::next(x, 1), std::next(dir, step * pos + offset1), E, Q, R);
      onestep_16(h2, h7, f2, *std::next(x, 2), std::next(dir, step * pos + offset2), E, Q, R);
      onestep_16(h3, h8, f3, *std::next(x, 3), std::next(dir, step * pos + offset3), E, Q, R);
      *std::next(hep, 2 * pos + 0) = h8;
      *std::next(hep, 2 * pos + 1) = E;
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
                           VECTORTYPE * Qm,
                           VECTORTYPE * Rm,
                           uint64_t ql,
                           VECTORTYPE * F0,
                           uint64_t * dir_long,
                           VECTORTYPE * H0,
                           VECTORTYPE * Mm,
                           VECTORTYPE * MQ,
                           VECTORTYPE * MR,
                           VECTORTYPE * MQ0) -> void
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
  for(auto pos = 0LL; pos < ql_signed; ++pos)
    {
      h4 = *std::next(hep, 2 * pos + 0);
      E  = *std::next(hep, 2 * pos + 1);
      VECTORTYPE * x = *std::next(qp, pos + 0);

      /* mask h4 and E */
      h4 = v_sub16(h4, *Mm);
      E  = v_sub16(E,  *Mm);

      /* init h4 and E */
      h4 = v_add16(h4, *MQ);
      E  = v_add16(E,  *MQ);
      E  = v_add16(E,  *MQ0);

      /* update MQ */
      *MQ = v_add16(*MQ,  *MR);

      onestep_16(h0, h5, f0, *std::next(x, 0), std::next(dir, step * pos + offset0), E, Q, R);
      onestep_16(h1, h6, f1, *std::next(x, 1), std::next(dir, step * pos + offset1), E, Q, R);
      onestep_16(h2, h7, f2, *std::next(x, 2), std::next(dir, step * pos + offset2), E, Q, R);
      onestep_16(h3, h8, f3, *std::next(x, 3), std::next(dir, step * pos + offset3), E, Q, R);
      *std::next(hep, 2 * pos + 0) = h8;
      *std::next(hep, 2 * pos + 1) = E;

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


auto search16(std::vector<WORD *> & q_start,
              WORD gap_open_penalty,
              WORD gap_extend_penalty,
              WORD * score_matrix,
              std::vector<WORD> & dprofile,
              WORD * hearray,
              uint64_t sequences,
              uint64_t * seqnos,
              uint64_t * scores,
              uint64_t * diffs,
              uint64_t * alignmentlengths,
              uint64_t qlen,
              std::vector<uint64_t> & dirbuffer) -> void
{
  static constexpr auto uint16_max = std::numeric_limits<uint16_t>::max();
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
  // an array of BYTES (or WORDS?)
  std::array<BYTE, channels * sizeof(VECTORTYPE) / sizeof(BYTE)> dseq {{}};

  uint64_t next_id {0};
  uint64_t done {0};

#ifdef __aarch64__
  const VECTORTYPE T0 = { uint16_max, 0, 0, 0, 0, 0, 0, 0 };
#elif defined __x86_64__
  const auto T0 = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, -1);
#elif defined __PPC__
  static constexpr auto unsigned_short_max = std::numeric_limits<unsigned short>::max();
  const VECTORTYPE T0 = { unsigned_short_max, 0, 0, 0, 0, 0, 0, 0 };
#endif

  assert((gap_open_penalty + gap_extend_penalty) <= std::numeric_limits<short>::max());
  assert(gap_extend_penalty <= std::numeric_limits<short>::max());
  auto Q = v_dup16(static_cast<short>(gap_open_penalty + gap_extend_penalty));
  auto R = v_dup16(static_cast<short>(gap_extend_penalty));

  done = 0;

  auto * hep = reinterpret_cast<VECTORTYPE *>(hearray);
  auto * * qp = reinterpret_cast<VECTORTYPE * *>(q_start.data());

  auto F0 = v_zero16();
  auto H0 = v_zero16();

  bool easy {false};

  uint64_t * dir = dirbuffer.data();

  while(true)
    {

      if (easy)
        {
          // fill all channels

          for(auto channel = 0U; channel < channels; ++channel)
            {
              for(auto j = 0U; j < cdepth; ++j)
                {
                  if (d_pos[channel] < d_length[channel]) {
                    dseq[channels * j + channel]
                      = 1 + nt_extract(d_address[channel], d_pos[channel]);
                    ++d_pos[channel];
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
              dprofile_shuffle16(dprofile.data(), score_matrix, dseq.data());
            }
          else
#endif
#endif
            {
              dprofile_fill16(dprofile.data(), score_matrix, dseq.data());
            }

#ifdef __x86_64__
#ifdef __SSE4_1__
          if (sse41_present != 0)
            {
              align_cells_regular_16_sse41(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
            }
          else
#endif
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

          M = v_zero16();
          T = T0;
          for(auto channel = 0U; channel < channels; ++channel)
            {
              if (d_pos[channel] < d_length[channel])
                {
                  // this channel has more sequence

                  for(auto j = 0U; j < cdepth; ++j)
                    {
                      if (d_pos[channel] < d_length[channel]) {
                        dseq[channels * j + channel]
                          = 1 + nt_extract(d_address[channel], d_pos[channel]);
                        ++d_pos[channel];
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

                  M = v_xor16(M, T);

                  const int64_t cand_id = seq_id[channel];

                  if (cand_id >= 0)
                    {
                      // save score

                      char * dbseq = d_address[channel];
                      const uint64_t dbseqlen = d_length[channel];
                      const uint64_t z = (dbseqlen + 3) % 4;
                      assert(z * channels + channel <= max_ptrdiff);
                      const uint64_t score
                        = *std::next(reinterpret_cast<WORD*>(S), static_cast<std::ptrdiff_t>(z * channels + channel));
                      *std::next(scores, cand_id) = score;

                      uint64_t diff {0};

                      if (score < uint16_max)
                        {
                          const uint64_t offset = d_offset[channel];
                          diff = backtrack<n_bits>(query.seq, dbseq, qlen, dbseqlen,
                                                   dirbuffer,
                                                   offset,
                                                   channel,
                                                   std::next(alignmentlengths, cand_id),
                                                   q_start.size());
                        }
                      else
                        {
                          diff = uint16_max;
                        }

                      *std::next(diffs, cand_id) = diff;

                      ++done;
                    }

                  if (next_id < sequences)
                    {
                      assert(next_id <= std::numeric_limits<int64_t>::max());
                      assert(next_id <= max_ptrdiff);
                      // get next sequence
                      seq_id[channel] = static_cast<int64_t>(next_id);
                      const uint64_t seqno = *std::next(seqnos, static_cast<std::ptrdiff_t>(next_id));
                      char * address {nullptr};
                      unsigned int length {0};

                      db_getsequenceandlength(seqno, address, length);

                      d_address[channel] = address;
                      d_length[channel] = length;

                      d_pos[channel] = 0;
                      d_offset[channel] = static_cast<uint64_t>(dir - dirbuffer.data());
                      ++next_id;

                      assert((2U * gap_open_penalty + 2U * gap_extend_penalty) <= std::numeric_limits<WORD>::max());
                      *std::next(reinterpret_cast<WORD *>(&H0), channel) = 0;
                      *std::next(reinterpret_cast<WORD *>(&F0), channel) = static_cast<WORD>(2U * gap_open_penalty + 2U * gap_extend_penalty);

                      // fill channel
                      for(auto j = 0U; j < cdepth; ++j)
                        {
                          if (d_pos[channel] < d_length[channel]) {
                            dseq[channels * j + channel] = 1 + nt_extract(d_address[channel], d_pos[channel]);
                            ++d_pos[channel];
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
                      for(auto j = 0U; j < cdepth; ++j) {
                        dseq[channels * j + channel] = 0;
                      }
                    }
                }

              T = v_shift_left16(T);
            }

          if (done == sequences) {
            break;
          }

#ifdef __x86_64__
#ifdef __SSE3__
          if (ssse3_present != 0)
            {
              dprofile_shuffle16(dprofile.data(), score_matrix, dseq.data());
            }
          else
#endif
#endif
            {
              dprofile_fill16(dprofile.data(), score_matrix, dseq.data());
            }

          MQ = v_and16(M, Q);
          MR = v_and16(M, R);
          MQ0 = MQ;

#ifdef __x86_64__
#ifdef __SSE4_1__
          if (sse41_present != 0)
            {
              align_cells_masked_16_sse41(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
            }
          else
#endif
#endif
            {
              align_cells_masked_16(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
            }
        }

      F0 = v_add16(F0, R);
      F0 = v_add16(F0, R);
      F0 = v_add16(F0, R);
      H0 = v_sub16(F0, Q);
      F0 = v_add16(F0, R);

      assert(4 * std::distance(q_start.begin(), q_start.end()) <= max_ptrdiff);
      dir = std::next(dir, 4 * std::distance(q_start.begin(), q_start.end()));
      auto const distance = std::distance(dirbuffer.begin(), dirbuffer.end());
      if (dir >= std::next(dirbuffer.data(), distance)) {
        dir = std::prev(dir, distance);
      }
    }
}
