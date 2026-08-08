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
#include "utils/search_data.hpp"  // Search_data (pulls in Cpu_features)
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
#include <vector>


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
static_assert(channels == channels_at_16_bits, "Dseq_16 is sized for 8 channels");
constexpr unsigned int cdepth {4};
constexpr uint8_t n_bits {16};
using BYTE = unsigned char;
using WORD = uint16_t;

// The three buffers arrive as their own containers rather than as three
// same-family pointers in a row: their types now say which is which, and
// the '_a'/'_v' parameters are unpacked once here into the base pointers
// the body walks. Inlining .data() at each cast site instead cost two
// prologue instructions and 16 bytes of search8 when the same thing was
// tried on hearray (commit 80df510).
auto dprofile_fill16(std::vector<WORD> & dprofile_v,
                            Score_matrix_16 const & score_matrix_a,
                            Dseq_16 const & dseq_a) -> void
{
  auto * const dprofile_word = dprofile_v.data();
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

inline auto onestep_16(VECTORTYPE & H,
                       VECTORTYPE & N,
                       VECTORTYPE & F,
                       VECTORTYPE const V,
                       WORD * const DIR,
                       VECTORTYPE & E,
                       VECTORTYPE const QR,
                       VECTORTYPE const R) -> void
{
  H = v_add16(H, V);
  auto W = H;
  H = v_min16(H, F);
  DIR[0] = v_mask_eq16(W, H);  // subscript, not std::next: hot loop, see align_cells
  H = v_min16(H, E);
  DIR[1] = v_mask_eq16(H, E);
  N = H;
  H = v_add16(H, QR);
  F = v_add16(F, R);
  E = v_add16(E, R);
  F = v_min16(H, F);
  DIR[2] = v_mask_eq16(H, F);
  E = v_min16(H, E);
  DIR[3] = v_mask_eq16(H, E);
}


// The masking payload: 'mask' selects the channels whose sequence just
// ended, 'mq' is the running gap-open accumulator seeded by the caller,
// 'mr' its per-iteration increment, and 'mq0' the value 'mq' held on
// entry. A plain struct rather than a template on VECTORTYPE: see
// utils/mask_vectors.hpp for why the template form is not usable here.
struct Mask_vectors {
  VECTORTYPE mask;
  VECTORTYPE mq;
  VECTORTYPE mr;
  VECTORTYPE mq0;
};


// The masking step, selected by the type of the kernel's mask argument
// (see utils/mask_vectors.hpp). The No_mask overload is empty, so the
// regular kernel's loop body contains nothing at this point.
inline auto apply_mask(VECTORTYPE &, VECTORTYPE &, No_mask const &) -> void
{
}

inline auto apply_mask(VECTORTYPE & h4, VECTORTYPE & E,
                       Mask_vectors & masks) -> void
{
  /* mask h4 and E */
  h4 = v_sub16(h4, masks.mask);
  E  = v_sub16(E,  masks.mask);

  /* init h4 and E */
  h4 = v_add16(h4, masks.mq);
  E  = v_add16(E,  masks.mq);
  E  = v_add16(E,  masks.mq0);

  /* update MQ */
  masks.mq = v_add16(masks.mq,  masks.mr);
}


// One block of cells, shared by the regular and masked kernels. The
// masked variant differs only by a per-iteration adjustment of h4 and E;
// which flavour this is comes from the type of 'masks', so the regular
// instantiation drops that adjustment entirely and is handed no masking
// data at all (see utils/mask_vectors.hpp).
template <typename Masks>
auto align_cells_16(VECTORTYPE * const Sm,
                    VECTORTYPE * const hep,
                    VECTORTYPE ** const qp,
                    VECTORTYPE const & Qm,
                    VECTORTYPE const & Rm,
                    uint64_t const ql,
                    VECTORTYPE const & F0,
                    uint64_t * const dir_long,
                    VECTORTYPE const & H0,
                    Masks & masks) -> void
{
  static constexpr auto step = 16;
  static constexpr auto offset0 = 0;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * dir = reinterpret_cast<WORD *>(dir_long);

  const auto Q = Qm;
  const auto R = Rm;

  auto f0 = F0;
  auto f1 = v_add16(f0, R);
  auto f2 = v_add16(f1, R);
  auto f3 = v_add16(f2, R);

  auto h0 = H0;
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
  // Performance: subscript / &dir[...] rather than std::next() in this hot loop
  // (same regression as commit 8c6925f on the SSE4.1 kernel). Stays clang-tidy
  // clean: pos is signed and operator[] is not pointer arithmetic.
  for (auto pos = 0LL; pos < ql_signed; ++pos)
    {
      VECTORTYPE const * x = qp[pos];
      h4 = hep[(2 * pos) + 0];
      E  = hep[(2 * pos) + 1];

      apply_mask(h4, E, masks);

      onestep_16(h0, h5, f0, x[0], &dir[(step * pos) + offset0], E, Q, R);
      onestep_16(h1, h6, f1, x[1], &dir[(step * pos) + offset1], E, Q, R);
      onestep_16(h2, h7, f2, x[2], &dir[(step * pos) + offset2], E, Q, R);
      onestep_16(h3, h8, f3, x[3], &dir[(step * pos) + offset3], E, Q, R);
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
  align_cells_16(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, no_mask);
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
  align_cells_16(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0, masks);
  *MQ = masks.mq;
}


namespace {

// Store the final score for the sequence that just ended in 'channel'
// and, when the score fits in a WORD, recover its number of differences
// by backtracking the alignment.
auto save_score_16(int64_t const cand_id,
                          unsigned int const channel,
                          VECTORTYPE const * const score_vectors,
                          std::array<Sequence, channels> const & d_sequence,
                          std::array<uint64_t, channels> const & d_offset,
                          Sequence const & query,
                          std::vector<uint64_t> const & dirbuffer,
                          uint64_t const q_start_size,
                          Span<uint64_t> const scores,
                          Span<uint64_t> const diffs,
                          uint64_t & done) -> void
{
  static constexpr auto uint16_max = std::numeric_limits<uint16_t>::max();

  // save score

  auto const & dbseq = d_sequence[channel];
  const uint64_t dbseqlen = dbseq.length;
  const uint64_t z = (dbseqlen + 3) % 4;
  assert(z * channels + channel <= max_ptrdiff);
  const uint64_t score
    = *std::next(reinterpret_cast<WORD const *>(score_vectors), static_cast<std::ptrdiff_t>((z * channels) + channel));
  assert(cand_id >= 0);
  auto const candidate = static_cast<std::size_t>(cand_id);
  scores[candidate] = score;

  uint64_t diff {0};

  if (score < uint16_max)
    {
      const uint64_t offset = d_offset[channel];
      diff = backtrack<n_bits>(query, dbseq,
                               make_view(dirbuffer),
                               offset,
                               channel,
                               q_start_size);
    }
  else
    {
      diff = uint16_max;
    }

  diffs[candidate] = diff;

  ++done;
}


// Write one 16-bit lane of a vector register.
//
// std::memcpy rather than a store through reinterpret_cast<WORD *>(&vec):
// a narrow store into an object whose declared type is VECTORTYPE is not
// something -fstrict-aliasing has to honour, so GCC is free to keep a
// stale copy of the vector in a register across it. This is not
// theoretical: with GCC 13.3 at -O3, 'swarm -d 4 -g 60' (a gap-open
// penalty high enough to select 16-bit mode at a low d, see
// set_bit_mode) produced clusters that disagreed with the -O0 build,
// and -fno-strict-aliasing alone restored them. memcpy aliases
// everything, so the lane write is always observed.
//
// Cold path: runs once per channel swap, never inside the kernel loop
// (measured free on 'd = 16', 18SV9-derived input).
auto set_lane_16(VECTORTYPE & vec, unsigned int const channel, WORD const value) -> void
{
  std::array<WORD, channels> lanes {{}};
  std::memcpy(lanes.data(), &vec, sizeof(vec));
  lanes[channel] = value;
  std::memcpy(&vec, lanes.data(), sizeof(vec));
}


// Attach the next database sequence to 'channel': record its address and
// length, reset the per-channel cursors, seed the H0/F0 lanes, and prime
// the first block. Returns whether the channel already reached the end of
// its (short) sequence, i.e. the next block is no longer "easy".
template <std::size_t capacity>
auto load_next_sequence_16(unsigned int const channel,
                                  Data const & data,
                                  View<uint64_t> const seqnos,
                                  uint64_t & next_id,
                                  uint64_t const * const dirbuffer_begin,
                                  uint64_t const * const dir,
                                  WORD const gap_open_penalty,
                                  WORD const gap_extend_penalty,
                                  VECTORTYPE & H0,
                                  VECTORTYPE & F0,
                                  std::array<unsigned char, capacity> & dseq,
                                  std::array<int64_t, channels> & seq_id,
                                  std::array<Sequence, channels> & d_sequence,
                                  std::array<uint64_t, channels> & d_pos,
                                  std::array<uint64_t, channels> & d_offset) -> bool
{
  assert(next_id <= std::numeric_limits<int64_t>::max());
  // get next sequence
  seq_id[channel] = static_cast<int64_t>(next_id);
  const uint64_t seqno = seqnos[next_id];
  auto const sequence = data.sequence_view(seqno);

  d_sequence[channel] = sequence;

  d_pos[channel] = 0;
  d_offset[channel] = static_cast<uint64_t>(dir - dirbuffer_begin);
  ++next_id;

  assert(((2U * gap_open_penalty) + (2U * gap_extend_penalty)) <= std::numeric_limits<WORD>::max());
  set_lane_16(H0, channel, 0);
  set_lane_16(F0, channel, static_cast<WORD>((2U * gap_open_penalty) + (2U * gap_extend_penalty)));

  // fill channel
  return fill_channel<channels, cdepth>(dseq, channel, d_sequence, d_pos);
}

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
  auto & dprofile = search_data.dprofile_w_v;
  auto * const hearray = search_data.hearray_v.data();  // He_block *
  auto const sequences = seqnos.size();
  auto const qlen = static_cast<uint64_t>(query.length);
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
  std::array<Sequence, channels> d_sequence {{}};
  std::array<int64_t, channels> seq_id {{}};
  seq_id.fill(-1);

  // refactoring fail: std::array -> warning: ignoring attributes on
  // template argument ‘VECTORTYPE’ {aka ‘__m128i’}
  VECTORTYPE S[4];

  // make an array of size VECTORTYPE * channels, but interpret as
  // an array of BYTES (or WORDS?)
  Dseq_16 dseq {{}};

  uint64_t next_id {0};
  uint64_t done {0};

  const auto T0 = make_T0_16();

  assert((gap_open_penalty + gap_extend_penalty) <= std::numeric_limits<short>::max());
  assert(gap_extend_penalty <= std::numeric_limits<short>::max());
  auto Q = v_dup16(static_cast<short>(gap_open_penalty + gap_extend_penalty));
  auto R = v_dup16(static_cast<short>(gap_extend_penalty));

  // one cast, from the over-aligned He_block straight to the vector type.
  // hearray used to be a WORD * that nothing else read, and that step
  // discarded the alignment these loads need.
  auto * hep = reinterpret_cast<VECTORTYPE *>(hearray);
  auto * * qp = reinterpret_cast<VECTORTYPE * *>(q_start.data());

  auto F0 = v_zero16();
  auto H0 = v_zero16();

  bool easy {false};

  uint64_t * dir = dirbuffer.data();

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

                  const int64_t cand_id = seq_id[channel];

                  if (cand_id >= 0)
                    {
                      save_score_16(cand_id, channel, S,
                                    d_sequence, d_offset,
                                    query, dirbuffer, q_start.size(),
                                    scores, diffs, done);
                    }

                  if (next_id < sequences)
                    {
                      if (load_next_sequence_16(channel, data, seqnos, next_id,
                                                dirbuffer.data(), dir,
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
      assert(dirbuffer.size() <= max_ptrdiff);
      if (dir >= std::next(dirbuffer.data(), static_cast<std::ptrdiff_t>(dirbuffer.size()))) {
        dir = std::prev(dir, static_cast<std::ptrdiff_t>(dirbuffer.size()));
      }
    }
}
