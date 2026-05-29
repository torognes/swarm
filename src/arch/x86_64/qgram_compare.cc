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

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#endif

#include "../../utils/qgram_compare.h"
#include "../../utils/cpu_features.h"  // Cpu_features
#include "../../utils/qgram_array.h"   // qgramvectorbytes
#include "popcnt.h"                    // compareqgramvectors_popcnt
#include <cassert>
#include <cstdint>  // uint64_t
#include <iterator>  // std::next


// C++20 refactoring: this SSE2-without-POPCNT path exists only for pre-Nehalem
// CPUs (before 2008). With std::popcount the dispatch (popcnt vs. SSE2 here)
// collapses into a single portable loop that the compiler lowers to the best
// available instruction for the target -march.

namespace {

auto v_add64(__m128i lhs, __m128i rhs) -> __m128i {
  // add 64-bit integers packed in lhs and rhs (SSE2)
  return _mm_add_epi64(lhs, rhs);
}

auto popcount_128(__m128i input_vector) -> uint64_t
{
  static constexpr unsigned char char1 {0x55};  // '0101 0101'
  static constexpr unsigned char char2 {0x33};  // '0011 0011'
  static constexpr unsigned char char4 {0x0f};  // '0000 1111'
  static constexpr auto shift_by_1 = 1;
  static constexpr auto shift_by_2 = 2;
  static constexpr auto shift_by_4 = 4;
  static constexpr auto shift_by_8 = 8;

  const auto mask1 = _mm_set_epi8(char1, char1, char1, char1, char1, char1, char1, char1,
                                  char1, char1, char1, char1, char1, char1, char1, char1);

  const auto mask2 = _mm_set_epi8(char2, char2, char2, char2, char2, char2, char2, char2,
                                  char2, char2, char2, char2, char2, char2, char2, char2);

  const auto mask4 = _mm_set_epi8(char4, char4, char4, char4, char4, char4, char4, char4,
                                  char4, char4, char4, char4, char4, char4, char4, char4);

  const auto zero = _mm_setzero_si128();

  /* add together 2 bits: 0+1, 2+3, 3+4, ... 126+127 */

  const auto vector_a = _mm_srli_epi64(input_vector, shift_by_1);
  const auto vector_b = _mm_and_si128(input_vector, mask1);
  const auto vector_c = _mm_and_si128(vector_a, mask1);
  const auto vector_d = v_add64(vector_b, vector_c);

  /* add together 4 bits: (0+1)+(2+3), ... (124+125)+(126+127) */

  const auto vector_e = _mm_srli_epi64(vector_d, shift_by_2);
  const auto vector_f = _mm_and_si128(vector_d, mask2);
  const auto vector_g = _mm_and_si128(vector_e, mask2);
  const auto vector_h = v_add64(vector_f, vector_g);

  /* add together 8 bits: (0..3)+(4..7), ... (120..123)+(124..127) */

  const auto vector_i = _mm_srli_epi64(vector_h, shift_by_4);
  const auto vector_j = v_add64(vector_h, vector_i);
  const auto vector_k = _mm_and_si128(vector_j, mask4);

  /* add together 8 bytes: (0..63) and (64..127) */

  const auto vector_l = _mm_sad_epu8(vector_k, zero);

  /* add together 64-bit values into final 128 bit value */

  const auto vector_m = _mm_srli_si128(vector_l, shift_by_8);
  const auto vector_n = v_add64(vector_m, vector_l);

  /* return low 64 bits: return value is always in range 0 to 128 */

  return reinterpret_cast<uint64_t>(_mm_movepi64_pi64(vector_n));
}


auto compareqgramvectors_128(unsigned char const * lhs, unsigned char const * rhs) -> uint64_t
{
  /* Count number of different bits */
  /* Uses SSE2 but not POPCNT instruction */
  assert(qgramvectorbytes % 16 == 0); // input MUST be 16-byte aligned

  static constexpr auto n_vector_lengths = qgramvectorbytes / sizeof(__m128i);  // 8
  auto const * lhs_ptr = reinterpret_cast<__m128i const *>(lhs);
  auto const * rhs_ptr = reinterpret_cast<__m128i const *>(rhs);
  uint64_t count {0};

  for (auto i = 0ULL; i < n_vector_lengths; ++i) {
    count += popcount_128(_mm_xor_si128(*lhs_ptr, *rhs_ptr));
    lhs_ptr = std::next(lhs_ptr);
    rhs_ptr = std::next(rhs_ptr);
  }

  return count;
}

}  // namespace


auto compareqgramvectors(unsigned char const * lhs, unsigned char const * rhs,
                         Cpu_features const & cpu_features) -> uint64_t
{
  if (cpu_features.popcnt) {
    return compareqgramvectors_popcnt(lhs, rhs);
  }
  return compareqgramvectors_128(lhs, rhs);
}
