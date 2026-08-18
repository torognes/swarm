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

#ifndef SWARM_UTILS_BACKTRACK_H
#define SWARM_UTILS_BACKTRACK_H

#include "../db.hpp"  // Sequence
#include "ceil_divide.hpp"  // ceil_divide
#include "nt_codec.hpp"
#include "view.hpp"  // View
#include <cassert>
#include <cstdint>  // uint64_t, uint8_t

#ifndef NDEBUG
#include <limits>
#endif


// The two search widths served by this header, named once here: the
// specialization below, backtrack()'s static_assert and search8.cpp /
// search16.cpp all spell the same fact.
constexpr uint8_t bits8 {8};
constexpr uint8_t bits16 {16};


// default template (16 bits)
template <uint8_t n_bits>
constexpr auto compute_mask(uint64_t const channel,
                            unsigned int const offset) -> uint64_t {
  return (3ULL << ((2 * channel) + offset));
}


// specialization (8 bits)
//
// Here rather than in search8.cpp, which used to declare it. An explicit
// specialization only applies where it has been declared, so a second
// translation unit instantiating backtrack<bits8> without seeing it would
// silently get the primary template instead -- and the two masks differ,
// so that miscomputes rather than failing to link. Declaring it beside
// the primary means every user of the header sees both.
//
// Definable in a header because it is constexpr, hence implicitly inline:
// an explicit specialization is not inline on its own, and without that
// specifier this definition would collide in every TU that included it.
template <>
constexpr auto compute_mask<bits8>(uint64_t const channel,
                                   unsigned int const offset) -> uint64_t {
  return (1ULL << (channel + offset));
}


// Both widths are constexpr, so which one a given n_bits selects is
// checkable at compile time -- and that is the thing the move above is
// for. These are what a translation unit that could not see the
// specialization would fail: it would resolve compute_mask<bits8> to the
// primary template and produce the 16-bit masks. (Same idea as the
// static_assert block in ceil_divide.hpp.)
static_assert(compute_mask<bits16>(0, 0) == 3ULL, "16 bits: two adjacent bits at channel 0");
static_assert(compute_mask<bits16>(1, 0) == 12ULL, "16 bits: a channel strides by two bits");
static_assert(compute_mask<bits8>(0, 0) == 1ULL, "8 bits: one bit at channel 0");
static_assert(compute_mask<bits8>(1, 0) == 2ULL, "8 bits: a channel strides by one bit");
static_assert(compute_mask<bits8>(1, 16) == (1ULL << 17), "8 bits: offset selects the mask field");

enum struct Alignment: unsigned char { Insertion, Deletion, Match };


// qseq and dseq each carry their own nucleotide count, so there is no
// way to pair one sequence's data with the other's length: the two
// parameters used to be four, adjacent and same-typed.
//
// Returns the number of differences in the optimal alignment. The
// alignment's length is computed on the way (it is what the difference
// count is derived from) but not reported: no caller reads it.
template <uint8_t n_bits>
auto backtrack(Sequence const & qseq,
               Sequence const & dseq,
               View<uint64_t> const dirbuffer,
               uint64_t const offset,
               uint64_t const channel,
               uint64_t const longestdbsequence) -> uint64_t {
  static_assert(n_bits == bits8 or n_bits == bits16, "n_bits must be 8 or 16");
  static constexpr auto offset0 = 0U;
  static constexpr auto offset1 = offset0 + 16;
  static constexpr auto offset2 = offset1 + 16;
  static constexpr auto offset3 = offset2 + 16;
  // refactoring C++17: if constexpr
  auto const maskup      = compute_mask<n_bits>(channel, offset0);
  auto const maskleft    = compute_mask<n_bits>(channel, offset1);
  auto const maskextup   = compute_mask<n_bits>(channel, offset2);
  auto const maskextleft = compute_mask<n_bits>(channel, offset3);

  // nucleotide counts; the packed words below are read by nucleotide_at
  auto const qlen = static_cast<uint64_t>(qseq.length);
  auto const dlen = static_cast<uint64_t>(dseq.length);
  assert(qlen <= std::numeric_limits<int64_t>::max());
  assert(dlen <= std::numeric_limits<int64_t>::max());
  auto column = static_cast<int64_t>(qlen) - 1;
  auto row = static_cast<int64_t>(dlen) - 1;
  uint64_t aligned {0};
  uint64_t matches {0};
  auto operation = Alignment::Match;  // no extension in progress yet

  // The aligner computes four rows per block, so the direction buffer
  // holds each block as 'longestdbsequence' columns of four sub-rows;
  // the same 4 spells the block's height, the column stride and the
  // sub-row within it.
  static constexpr uint64_t rows_per_block {4};

  // That buffer is a ring: the aligner advances its write cursor by one
  // block per iteration and wraps it by subtracting the size once
  // (search8.cpp:876-878, search16.cpp:624-626). The read below wraps the
  // same way, which is exact only while every index stays under twice the
  // size:
  //  - 'offset' is a cursor position, recorded before that cursor's own
  //    wrap, so it is inside the buffer;
  //  - the largest cell this loop can address is the one for the last row
  //    and the last column, which is below
  //    rows_per_block * longestdbsequence * ceil(dlen / rows_per_block) --
  //    and scanner.cpp:56 allocates exactly that with both lengths at
  //    their maximum, since a query and a database sequence are both
  //    database sequences.
  // Their sum is therefore below twice the size. Wrapping with '%' instead
  // would be a 64-bit hardware division on every iteration of this loop,
  // which is the most expensive thing in it by an order of magnitude.
  auto const ring_size = dirbuffer.size();
  assert(offset < ring_size);
  assert(qlen <= longestdbsequence);
  assert(dlen <= longestdbsequence);
  assert(rows_per_block * longestdbsequence * ceil_divide(dlen, rows_per_block)
         <= ring_size);

  while ((column >= 0) and (row >= 0)) {
      ++aligned;

      auto const row_index = static_cast<uint64_t>(row);
      auto const cell
        = (longestdbsequence * rows_per_block * (row_index / rows_per_block))
        + (rows_per_block * static_cast<uint64_t>(column))
        + (row_index % rows_per_block);
      auto index = offset + cell;
      if (index >= ring_size) { index -= ring_size; }
      auto const direction = dirbuffer[index];

      if ((operation == Alignment::Insertion) and ((direction & maskextleft) == 0U)) {
        --row;
      }
      else if ((operation == Alignment::Deletion) and ((direction & maskextup) == 0U)) {
        --column;
      }
      else if ((direction & maskleft) != 0U) {
        --row;
        operation = Alignment::Insertion;
      }
      else if ((direction & maskup) == 0U) {
        --column;
        operation = Alignment::Deletion;
      }
      else
        {
          if (nucleotide_at(qseq, static_cast<uint64_t>(column)) ==
              nucleotide_at(dseq, static_cast<uint64_t>(row))) {
            ++matches;
          }
          --column;
          --row;
          operation = Alignment::Match;
        }
    }

  // The loop above stopped as soon as one of the two ran out, so at most
  // one of these remainders is non-negative; each position still left is
  // one more aligned column.
  if (column >= 0) { aligned += static_cast<uint64_t>(column) + 1; }
  if (row >= 0)    { aligned += static_cast<uint64_t>(row) + 1; }

  return aligned - matches;
}

#endif  // SWARM_UTILS_BACKTRACK_H
