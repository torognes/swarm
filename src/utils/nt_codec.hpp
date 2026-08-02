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

#ifndef SWARM_UTILS_NT_CODEC_H
#define SWARM_UTILS_NT_CODEC_H

#include <cstdint>  // uint64_t


// Nucleotides are packed four per byte, as 2-bit fields, lowest field
// first. Locating the byte and decoding within it are two separate steps:
//
//   nt_byte_index(position)      which packed byte holds a nucleotide
//   nt_extract(byte, position)   decode that nucleotide out of the byte
//
// nt_extract() therefore touches no memory and needs no length: it is a
// function of one byte and one position. The caller does the lookup, which
// is where the bound is known -- and where a loop walking a sequence in
// order can load one byte and decode four nucleotides from it instead of
// re-loading the same byte four times. Whether that last trade pays is a
// question for the caller and not for this header: fill_channel()
// (dseq_fill.hpp) tried it and measured slower.
//
// The same field order describes the buffer at word granularity, which
// matters more than it looks: nucleotide 'position' occupies bits
// [2 * position, 2 * position + 2) of the packed bytes read as a
// little-endian bit array, so the byte view above and a 64-bit word view
// agree with no conversion between them. That is what lets
// seq_identical() (variants.cpp) compare 32 nucleotides with a single
// xor, and nt_bytelength() below rounds every sequence up to a whole
// multiple of 8 bytes, so such a word read always stays inside it.
//
// Do not flip the field order to make a packed byte read left to right.
// Reversing it inside the byte alone breaks the agreement between those
// two views, and the order is baked into the packer (db.cpp), nt_set()
// and the window helpers (variants.cpp), and the precomputed byte-rate
// table in zobrist.cpp: a site left out of step changes hashes, and
// therefore clusters, without failing.
//
// Callers holding a Sequence should prefer nucleotide_at() (db.hpp), which
// does both steps and checks the position against the sequence's
// nucleotide count.
//
// Both live here rather than in nt_codec.cpp so that they inline without
// relying on link-time optimisation.

inline auto nt_byte_index(uint64_t const position) -> uint64_t {
  // 4 nucleotides stored per byte, so nucleotide 34 lives in byte 8
  static constexpr auto divide_by_4 = 2U;
  return position >> divide_by_4;
}


inline auto nt_extract(char const compressed_byte, uint64_t const position) -> unsigned char {
  // Extract a given position from a compressed byte
  //
  // example: extract nucleotide at position 34
  //  - (note: coordinates are zero-based),
  //  - (note: 4 nucleotides stored per byte),
  //  - 34 & mask_upper_bits -> 2
  //    (nucleotide is stored in the pair of bits at position 2 in the
  //    compressed byte),
  //  - left-shift compressed byte 2 times (equivalent to dividing by 4),
  //  - the pair of bits we are looking for is now at the start of the byte,
  //  - mask upper bits to keep only the encoded nucleotide (-> 0, 1, 2, or 3)
  //
  static constexpr auto max_nt_per_byte = 4U; // 4 nt fit in 8 bits
  static constexpr auto keep_first_two_bits = max_nt_per_byte - 1; // 0000'0011 (mask all upper bits)
  auto const target_pair_of_bits = position & keep_first_two_bits;  // same as pos & 4 (remainder): 0, 1, 2, or 3
  auto const divider = target_pair_of_bits << 1U;  // left-shift by 0, 2, 4, or 6 (same as dividing by 0, 4, 16, or 64)

  // outputs four possible values: 0, 1, 2 or 3
  // (cast to unsigned int avoids the integer promotion that would
  //  otherwise apply the right shift to a signed int and trip the
  //  hicpp-signed-bitwise check)
  return static_cast<unsigned char>(
      (static_cast<unsigned int>(static_cast<unsigned char>(compressed_byte)) >> divider)
      & keep_first_two_bits);
}


auto nt_bytelength(unsigned int len) -> unsigned int;

#endif
