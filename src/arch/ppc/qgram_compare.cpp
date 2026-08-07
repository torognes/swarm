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

#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#include "../../utils/qgram_compare.hpp"
#include "../../utils/cpu_features.hpp"  // Cpu_features
#include "../../utils/qgram_array.hpp"   // qgramvectorbytes
#include <cstdint>  // uint64_t


// C++20 refactoring: replace with a portable loop using std::popcount; on
// ppc64le the compiler emits vpopcntd for tight XOR + popcount loops.
auto compareqgramvectors(Qgram_vector const & lhs, Qgram_vector const & rhs,
                         Cpu_features const & cpu_features) -> uint64_t
{
  static_cast<void>(cpu_features);  // unused: AltiVec vpopcnt is always available on Power8+
  static_assert(qgramvectorbytes % sizeof(vector unsigned char) == 0,
                "qgram vector must be a whole number of 128-bit words");
  static_assert(alignof(Qgram_vector) >= alignof(vector unsigned char),
                "qgram vector must be aligned for the loads below");
  static constexpr auto n_vector_lengths = qgramvectorbytes / sizeof(vector unsigned char);  // 8
  // cast from &lhs, not lhs.data(): data() hands back unsigned char const *,
  // which drops the alignas(16) these aligned loads rely on
  auto const * lhs_ptr = reinterpret_cast<vector unsigned char const *>(&lhs);
  auto const * rhs_ptr = reinterpret_cast<vector unsigned char const *>(&rhs);
  vector unsigned long long count_vector = { 0, 0 };

  for (auto i = 0ULL; i < n_vector_lengths; ++i) {
    count_vector += vec_vpopcnt(reinterpret_cast<vector unsigned long long>(vec_xor(*lhs_ptr, *rhs_ptr)));
    ++lhs_ptr;
    ++rhs_ptr;
  }

  return count_vector[0] + count_vector[1];
}
