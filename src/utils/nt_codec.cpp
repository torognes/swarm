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

#include "nt_codec.hpp"
#include <cassert>

#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
#include <limits>
constexpr auto uint_max = std::numeric_limits<unsigned int>::max();
#endif


// round-up operation, compiler cleverly elimates the multiplication (8 is a power of 2)
auto nt_bytelength(unsigned int const len) -> unsigned int {
  // Compute number of bytes used for compressed sequence of length len
  // (minimum result is 8 bytes)
  static constexpr auto max_nt_per_uint64 = 32U;  // 32 nt fit in 64 bits
  static constexpr auto divide_by_32 = 5U;  // (len + 31) % 32 (drop remainder)
  static constexpr auto bytes_per_uint64 = 8U;  // times 8 to get the number of bytes
  assert(len != 0);
  assert(len <= uint_max - (max_nt_per_uint64 - 1));
  return ((len + max_nt_per_uint64 - 1) >> divide_by_32) * bytes_per_uint64;
}
