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

// same aliases as intrinsics_to_functions.cpp; '__vector' is the
// standards-mode spelling of the AltiVec 'vector' keyword
using v_u8_t = __vector unsigned char;
using v_u64_t = __vector unsigned long long;


// C++20 refactoring: replace with a portable loop using std::popcount; on
// ppc64le the compiler emits vpopcntd for tight XOR + popcount loops.
auto compareqgramvectors(Qgram_vector const & lhs, Qgram_vector const & rhs,
                         Cpu_features const & cpu_features) -> uint64_t
{
  static_cast<void>(cpu_features);  // unused: AltiVec vpopcnt is always available on Power8+
  // false positive: cppcheck does not know the 16-byte __vector types and
  // evaluates sizeof(v_u8_t) as sizeof(unsigned char)
  // cppcheck-suppress moduloofone
  static_assert(qgramvectorbytes % sizeof(v_u8_t) == 0,
                "qgram vector must be a whole number of 128-bit words");
  static_assert(alignof(Qgram_vector) >= alignof(v_u8_t),
                "qgram vector must be aligned for the loads below");
  static constexpr auto n_vector_lengths = qgramvectorbytes / sizeof(v_u8_t);  // 8
  // cast from &lhs, not lhs.data(): data() hands back unsigned char const *,
  // which drops the alignas(16) these aligned loads rely on
  auto const * lhs_ptr = reinterpret_cast<v_u8_t const *>(&lhs);
  auto const * rhs_ptr = reinterpret_cast<v_u8_t const *>(&rhs);
  v_u64_t count_vector = { 0, 0 };

  for (auto i = 0ULL; i < n_vector_lengths; ++i) {
    count_vector += vec_vpopcnt(reinterpret_cast<v_u64_t>(vec_xor(*lhs_ptr, *rhs_ptr)));
    ++lhs_ptr;
    ++rhs_ptr;
  }

  return count_vector[0] + count_vector[1];
}
