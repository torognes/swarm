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

#include <cstdint>  // uint64_t
#include <iterator>  // std::next


inline auto nt_extract(char * seq, uint64_t pos) -> unsigned char
{
  // Extract compressed nucleotide in sequence seq at position pos
  // refactoring: replace with bitset operations?
  static constexpr auto max_nt_per_uint64 = 32U;  // 32 nt fit in 64 bits
  static constexpr auto drop_remainder = 5U;  // (len+31) % 32 (drop remainder)
  static constexpr auto max_range = 3U;
  const auto target_uint64 = static_cast<long long int>(pos >> drop_remainder);  // same as pos % 32
  auto const * seq_uint64 = reinterpret_cast<uint64_t *>(seq);
  const auto compressed_chunk = *std::next(seq_uint64, target_uint64);
  // max_nt_per_uint64 - 1 = 31 = ...OOO1 1111 (mask all upper bits)
  const auto target_nucleotide = (pos & (max_nt_per_uint64 - 1)) << 1U;  // 0, 2, 4,..., 60, 62
  // outputs four possible values: 0, 1, 2 or 3
  return (compressed_chunk >> target_nucleotide) & max_range;  // UBSAN: misaligned address for type 'uint64_t', which requires 8 byte alignment
}


inline auto nt_bytelength(unsigned int len) -> unsigned int
{
  // Compute number of bytes used for compressed sequence of length len
  // (minimum result is 8 bytes)
  static constexpr auto max_nt_per_uint64 = 32U;  // 32 nt fit in 64 bits
  static constexpr auto drop_remainder = 5U;  // (len + 31) % 32 (drop remainder)
  static constexpr auto bytes_per_uint64 = 8U;  // times 8 to get the number of bytes
  return ((len + max_nt_per_uint64 - 1) >> drop_remainder) * bytes_per_uint64;
}
