/*
    SWARM

    Copyright (C) 2012-2025 Torbjorn Rognes and Frederic Mahe

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

#include <cassert>
#include <cstdint>  // uint64_t
#include <iterator>  // std::next

#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
#include <limits>
constexpr auto uint_max = std::numeric_limits<unsigned int>::max();
#endif


auto nt_extract(char const * compressed_sequence, uint64_t const position) -> unsigned char {
  // Extract a given position from a compressed sequence
  //
  // example: extract nucleotide at position 34
  //  - (note: coordinates are zero-based),
  //  - (note: 4 nucleotides stored per byte),
  //  - 34 / 4 -> 8
  //    (nucleotide is stored in the byte at position 8 in the compressed
  //    sequence),
  //  - 34 & mask_upper_bits -> 2
  //    (nucleotide is stored in the pair of bits at position 2 in the
  //    compressed byte),
  //  - left-shift compressed byte 2 times (equivalent to dividing by 4),
  //  - the pair of bits we are looking for is now at the start of the byte,
  //  - mask upper bits to keep only the encoded nucleotide (-> 0, 1, 2, or 3)
  //
  static constexpr auto divide_by_4 = 2U;
  assert((position >> divide_by_4) <= std::numeric_limits<long long int>::max());
  auto const target_byte = static_cast<long long int>(position >> divide_by_4);  // same as int{pos / 4}
  auto const compressed_byte = static_cast<unsigned char>(*std::next(compressed_sequence, target_byte));
  static constexpr auto max_nt_per_byte = 4U; // 4 nt fit in 8 bits
  static constexpr auto keep_first_two_bits = max_nt_per_byte - 1; // 0000'0011 (mask all upper bits)
  auto const target_pair_of_bits = position & keep_first_two_bits;  // same as pos & 4 (remainder): 0, 1, 2, or 3
  auto const divider = target_pair_of_bits << 1U;  // left-shift by 0, 2, 4, or 6 (same as dividing by 0, 4, 16, or 64)

  // outputs four possible values: 0, 1, 2 or 3
  return (compressed_byte >> divider) & keep_first_two_bits;
}


auto nt_bytelength(const unsigned int len) -> unsigned int
{
  // Compute number of bytes used for compressed sequence of length len
  // (minimum result is 8 bytes)
  static constexpr auto max_nt_per_uint64 = 32U;  // 32 nt fit in 64 bits
  static constexpr auto divide_by_32 = 5U;  // (len + 31) % 32 (drop remainder)
  static constexpr auto bytes_per_uint64 = 8U;  // times 8 to get the number of bytes
  assert(len != 0);
  assert(len <= uint_max - (max_nt_per_uint64 - 1));
  return ((len + max_nt_per_uint64 - 1) >> divide_by_32) * bytes_per_uint64;
}
