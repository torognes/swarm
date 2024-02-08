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


inline auto nt_extract(char * seq, uint64_t pos) -> unsigned char
{
  // Extract compressed nucleotide in sequence seq at position pos
  static constexpr unsigned int max_nt_per_uint64 {32};  // 32 nt fit in 64 bits
  static constexpr unsigned int drop_remainder {5};  // (len+31) % 32 (drop remainder)
  static constexpr unsigned int max_range {3};
  // outputs four possible values: 0, 1, 2 or 3
  return (((reinterpret_cast<uint64_t*>(seq))[pos >> drop_remainder]) >> \
          ((pos & (max_nt_per_uint64 - 1)) << 1U)) & max_range;   // UBSAN: misaligned address for type 'uint64_t', which requires 8 byte alignment
}


inline auto nt_bytelength(unsigned int len) -> unsigned int
{
  // Compute number of bytes used for compressed sequence of length len
  // (minimum result is 8 bytes)
  static constexpr unsigned int max_nt_per_uint64 {32};  // 32 nt fit in 64 bits
  static constexpr unsigned int drop_remainder {5};  // (len + 31) % 32 (drop remainder)
  static constexpr unsigned int bytes_per_uint64 {8};  // times 8 to get the number of bytes
  return ((len + max_nt_per_uint64 - 1) >> drop_remainder) * bytes_per_uint64;
}
