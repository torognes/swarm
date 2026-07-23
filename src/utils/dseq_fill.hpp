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

#ifndef SWARM_UTILS_DSEQ_FILL_H
#define SWARM_UTILS_DSEQ_FILL_H

#include "nt_codec.hpp"  // nt_extract
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t


// Channel-filling helpers shared by search8 and search16. These operate
// only on the scalar 'dseq' staging buffer and the per-channel database
// cursors (d_address / d_pos / d_length); they do not touch any SIMD
// vector type, so they are identical for both search widths apart from
// the 'channels' (8 vs 16) and 'cdepth' constants, which are passed as
// template arguments.
//
// 'dseq' is passed as its std::array (its capacity is deduced) rather
// than as a raw pointer, so element access stays bounds-aware and free
// of pointer arithmetic. It is laid out as 'cdepth' blocks of 'channels'
// bytes each: the nucleotide at depth j for a given channel lives at
// index (channels * j) + channel.


// Fill the 'cdepth' depth-slots of a single channel from its database
// sequence, advancing d_pos. Slots past the end of the sequence are
// zero-filled. Returns true when the channel has reached the end of its
// sequence (d_pos == d_length), i.e. the block is no longer "easy".
template <unsigned int channels, unsigned int cdepth, std::size_t capacity>
inline auto fill_channel(std::array<unsigned char, capacity> & dseq,
                         unsigned int const channel,
                         std::array<char const *, channels> const & d_address,
                         std::array<uint64_t, channels> & d_pos,
                         std::array<uint64_t, channels> const & d_length) -> bool
{
  for (auto j = 0U; j < cdepth; ++j)
    {
      if (d_pos[channel] < d_length[channel]) {
        dseq[(channels * j) + channel]
          = 1 + nt_extract(d_address[channel], d_pos[channel]);
        ++d_pos[channel];
      }
      else {
        dseq[(channels * j) + channel] = 0;
      }
    }
  return d_pos[channel] == d_length[channel];
}


// Fill every channel (the "easy" block, where no channel switch is
// pending). Returns whether the next block is still easy: false as soon
// as any channel reaches the end of its sequence.
template <unsigned int channels, unsigned int cdepth, std::size_t capacity>
inline auto fill_all_channels(std::array<unsigned char, capacity> & dseq,
                              std::array<char const *, channels> const & d_address,
                              std::array<uint64_t, channels> & d_pos,
                              std::array<uint64_t, channels> const & d_length) -> bool
{
  bool easy {true};
  for (auto channel = 0U; channel < channels; ++channel)
    {
      if (fill_channel<channels, cdepth>(dseq, channel, d_address, d_pos, d_length)) {
        easy = false;
      }
    }
  return easy;
}


// Zero the 'cdepth' depth-slots of an empty channel (no more sequences
// to process in that channel).
template <unsigned int channels, unsigned int cdepth, std::size_t capacity>
inline auto clear_channel(std::array<unsigned char, capacity> & dseq,
                          unsigned int const channel) -> void
{
  for (auto j = 0U; j < cdepth; ++j) {
    dseq[(channels * j) + channel] = 0;
  }
}

#endif
