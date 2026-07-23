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

#include "../../utils/qgram_compare.hpp"
#include "../../utils/cpu_features.hpp"  // Cpu_features
#include "../../utils/qgram_array.hpp"   // qgramvectorbytes
#include <arm_neon.h>
#include <cstdint>  // uint64_t


// C++20 refactoring: replace with a portable loop using std::popcount; on
// aarch64 with -O3 the compiler auto-vectorizes to cnt + addv anyway.
auto compareqgramvectors(unsigned char const * lhs, unsigned char const * rhs,
                         Cpu_features const & cpu_features) -> uint64_t
{
  static_cast<void>(cpu_features);  // unused: NEON cnt is always available
  static constexpr auto n_vector_lengths = qgramvectorbytes / sizeof(uint8x16_t);  // 8
  auto const * lhs_ptr = reinterpret_cast<uint8x16_t const *>(lhs);
  auto const * rhs_ptr = reinterpret_cast<uint8x16_t const *>(rhs);
  uint64_t count {0};

  for (auto i = 0ULL; i < n_vector_lengths; ++i) {
    count += vaddvq_u8(vcntq_u8(veorq_u8(*lhs_ptr, *rhs_ptr)));
    ++lhs_ptr;
    ++rhs_ptr;
  }

  return count;
}
