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

#include "nt_codec.h"
#include <cstdint>  // uint64_t

// default template (16 bits)
template <uint8_t n_bits>
auto compute_mask(uint64_t const channel,
                  unsigned int const offset) -> uint64_t {
  return (3ULL << (2 * channel + offset));
}


// refactoring: could 'Unknown' be eliminated?
enum struct Alignment: unsigned char { Unknown, Insertion, Deletion, Match };


template <uint8_t n_bits>
auto backtrack(char * qseq,
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
  static constexpr uint8_t bits8 {8};
  static constexpr uint8_t bits16 {16};
  static_assert(n_bits == bits8 || n_bits == bits16, "n_bits must be 8 or 16");
  static constexpr auto offset0 {0U};
  static constexpr auto offset1 {offset0 + 16};
  static constexpr auto offset2 {offset1 + 16};
  static constexpr auto offset3 {offset2 + 16};
  const uint64_t maskup      = compute_mask<n_bits>(channel, offset0);
  const uint64_t maskleft    = compute_mask<n_bits>(channel, offset1);
  const uint64_t maskextup   = compute_mask<n_bits>(channel, offset2);
  const uint64_t maskextleft = compute_mask<n_bits>(channel, offset3);

  auto column = static_cast<int64_t>(qlen) - 1;
  auto row = static_cast<int64_t>(dlen) - 1;
  uint64_t aligned {0};
  uint64_t matches {0};
  Alignment operation = Alignment::Unknown;  // Insertion, Deletion or Match

  while ((column >= 0) and (row >= 0))
    {
      ++aligned;

      const uint64_t d
        = dirbuffer[(offset
                     + longestdbsequence * 4 * static_cast<uint64_t>(row / 4)
                     + 4 * static_cast<uint64_t>(column)
                     + (static_cast<uint64_t>(row) & 3U)
                     ) % dirbuffersize];  // refactoring: how to rename that variable?

      if ((operation == Alignment::Insertion) and ((d & maskextleft) == 0U))
        {
          --row;
        }
      else if ((operation == Alignment::Deletion) and ((d & maskextup) == 0U))
        {
          --column;
        }
      else if ((d & maskleft) != 0U)
        {
          --row;
          operation = Alignment::Insertion;
        }
      else if ((d & maskup) == 0U)
        {
          --column;
          operation = Alignment::Deletion;
        }
      else
        {
          if (nt_extract(qseq, static_cast<uint64_t>(column)) ==
              nt_extract(dseq, static_cast<uint64_t>(row))) {
            ++matches;
          }
          --column;
          --row;
          operation = Alignment::Match;
        }
    }

  while (column >= 0)
    {
      ++aligned;
      --column;
    }

  while (row >= 0)
    {
      ++aligned;
      --row;
    }

  * alignmentlengthp = aligned;
  return aligned - matches;
}
