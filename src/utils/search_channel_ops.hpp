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


#ifndef SWARM_UTILS_SEARCH_CHANNEL_OPS_H
#define SWARM_UTILS_SEARCH_CHANNEL_OPS_H


// What happens when a channel finishes a database sequence: store its
// score (and recover its difference count), then load the next sequence
// into that channel, seeding the two lane registers.
//
// These were three functions written twice, once per search width, and the
// two copies differed only in which width they named -- the lane type
// (BYTE against WORD), the saturation ceiling, the channel count and the
// backtrack width, each of which is a compile-time fact already recorded
// in utils/search_data.hpp. The copies had drifted: set_lane's comment
// carried the full account of a miscompilation on one side and a summary
// on the other, and load_next_sequence's assertions had moved around.
//
// Written against names the including translation unit provides rather
// than as a template on the lane type, for the same reason
// utils/align_cells_16.hpp is: VECTORTYPE cannot be a template argument
// without GCC dropping its attributes, may_alias among them (see
// utils/mask_vectors.hpp), and set_lane takes a VECTORTYPE &. Given that
// one name has to come from the includer, they all may, and the parameter
// lists stay readable.
//
// So this header includes nothing, and must be included last, inside the
// includer's anonymous namespace, with all of these in scope:
//
//   LANE           the lane type: BYTE at 8 bits, WORD at 16
//   score_ceiling  the value a saturating lane sticks at, for that LANE
//   VECTORTYPE     the channel register type, per architecture and width
//   channels       lanes per register, cdepth   depth slots per block
//   n_bits         the backtrack width, 8 or 16
//   max_ptrdiff    used by the assertions, debug builds only


// Store the final score for the sequence that just ended in 'channel'
// and, when the score fits in a LANE, recover its number of differences
// by backtracking the alignment.
auto save_score(int64_t const cand_id,
                unsigned int const channel,
                VECTORTYPE const * const score_vectors,
                std::array<Sequence, channels> const & d_sequence,
                std::array<uint64_t, channels> const & d_offset,
                Sequence const & query,
                View<uint64_t> const dirbuffer,
                uint64_t const q_start_size,
                Span<uint64_t> const scores,
                Span<uint64_t> const diffs,
                uint64_t & done) -> void
{
  // save score

  auto const & dbseq = d_sequence[channel];
  uint64_t const dbseqlen = dbseq.length;
  uint64_t const z = (dbseqlen + 3) % 4;
  assert(z * channels + channel <= max_ptrdiff);
  uint64_t const score
    = *std::next(reinterpret_cast<LANE const *>(score_vectors), static_cast<std::ptrdiff_t>((z * channels) + channel));
  assert(cand_id >= 0);
  auto const candidate = static_cast<std::size_t>(cand_id);
  scores[candidate] = score;

  uint64_t diff {0};

  if (score < score_ceiling)
    {
      uint64_t const offset = d_offset[channel];
      diff = backtrack<n_bits>(query, dbseq,
                               dirbuffer,
                               offset,
                               channel,
                               q_start_size);
    }
  else
    {
      diff = score_ceiling;
    }

  diffs[candidate] = diff;

  ++done;
}


// Write one lane of a vector register.
//
// std::memcpy rather than a store through reinterpret_cast<LANE *>(&vec):
// a narrow store into an object whose declared type is VECTORTYPE is not
// something -fstrict-aliasing has to honour, so GCC is free to keep a
// stale copy of the vector in a register across it and drop the write.
// This is not theoretical: with GCC 13.3 at -O3, 'swarm -d 4 -g 60' (a
// gap-open penalty high enough to select 16-bit mode at a low d, see
// set_bit_mode) produced clusters that disagreed with the -O0 build, and
// -fno-strict-aliasing alone restored them. memcpy aliases everything, so
// the lane write is always observed. It was the 16-bit width that was
// caught misbehaving; the 8-bit one never was, but the construct is the
// same one and now so is the code.
//
// Cold path: runs once per channel swap, never inside the kernel loop
// (measured free on 'd = 16', 18SV9-derived input).
auto set_lane(VECTORTYPE & vec, unsigned int const channel, LANE const value) -> void
{
  std::array<LANE, channels> lanes {{}};
  std::memcpy(lanes.data(), &vec, sizeof(vec));
  lanes[channel] = value;
  std::memcpy(&vec, lanes.data(), sizeof(vec));
}


// Hand 'channel' the next database sequence, seeding its H0 and F0 lanes
// and filling the first block. Returns whether the channel already reached
// the end of its (short) sequence, i.e. the next block is no longer "easy".
template <std::size_t capacity>
auto load_next_sequence(unsigned int const channel,
                        Data const & data,
                        View<uint64_t> const seqnos,
                        uint64_t & next_id,
                        Span<uint64_t> const dirbuffer,
                        uint64_t const * const dir,
                        LANE const gap_open_penalty,
                        LANE const gap_extend_penalty,
                        VECTORTYPE & H0,
                        VECTORTYPE & F0,
                        std::array<unsigned char, capacity> & dseq,
                        std::array<int64_t, channels> & seq_id,
                        std::array<Sequence, channels> & d_sequence,
                        std::array<uint64_t, channels> & d_pos,
                        std::array<uint64_t, channels> & d_offset) -> bool
{
  // get next sequence
  assert(next_id <= std::numeric_limits<int64_t>::max());
  seq_id[channel] = static_cast<int64_t>(next_id);
  uint64_t const seqno = seqnos[next_id];
  auto const sequence = data.sequence_view(seqno);
  d_sequence[channel] = sequence;
  d_pos[channel] = 0;
  // the direction-buffer cursor this channel's backtrack will start from,
  // as an offset rather than a pointer: save_score turns it into a subview
  d_offset[channel] = static_cast<uint64_t>(dir - dirbuffer.cbegin());
  ++next_id;

  set_lane(H0, channel, 0);
  assert((2U * gap_open_penalty) + (2U * gap_extend_penalty) <= std::numeric_limits<LANE>::max());
  set_lane(F0, channel, static_cast<LANE>((2U * gap_open_penalty) + (2U * gap_extend_penalty)));

  // fill channel
  return fill_channel<channels, cdepth>(dseq, channel, d_sequence, d_pos);
}

#endif  // SWARM_UTILS_SEARCH_CHANNEL_OPS_H
