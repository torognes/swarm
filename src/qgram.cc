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

#include "swarm.h"

#ifdef __aarch64__
#include <arm_neon.h>
#endif


#ifdef __PPC__
#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif
#endif


#ifdef __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#include "utils/x86_cpu_feature_popcnt.h"
#endif

#include "popcnt.h"  // refactoring: fence with #ifdef __POPCNT__?

#endif

#include "utils/qgram_array.h"
#include "utils/qgram_threadinfo.h"
#include "utils/nt_codec.h"
#include "utils/threads.h"
#include <cassert>
#include <cstdint>  // int64_t, uint64_t
#include <cstring>  // memset
#include <iterator>  // std::next
#include <limits>
#include <vector>


qgramvector_t * qgrams {nullptr};
static ThreadRunner * qgram_threads = nullptr;
static struct thread_info_s * thread_info_ptr;


auto findqgrams(unsigned char * seq, uint64_t seqlen,
                unsigned char * qgramvector) -> void
{
  /* set qgram bit vector by xoring occurrences of qgrams in sequence */

  static constexpr unsigned int max_range {7};

  std::memset(qgramvector, 0, qgramvectorbytes);

  uint64_t qgram {0};
  unsigned int position {0};

  while((position < qgramlength - 1) and (position < seqlen))
  {
    qgram = (qgram << 2U) | nt_extract(reinterpret_cast<char *>(seq), position);
    ++position;
  }

  while(position < seqlen)
  {
    qgram = (qgram << 2U) | nt_extract(reinterpret_cast<char *>(seq), position);
    assert((qgram & max_range) <= 7);
    qgramvector[(qgram >> 3U) & (qgramvectorbytes - 1)] ^= static_cast<unsigned char>(1U << (qgram & max_range));
    ++position;
  }
}

auto qgram_work_diff(thread_info_s * tip) -> void;
auto qgram_worker(int64_t nth_thread) -> void;
auto compareqgramvectors(unsigned char * lhs, unsigned char * rhs) -> uint64_t;

#ifdef __aarch64__

auto compareqgramvectors(unsigned char * lhs, unsigned char * rhs) -> uint64_t
{
  static constexpr auto n_vector_lengths = qgramvectorbytes / sizeof(uint8x16_t);  // 8
  auto * lhs_ptr = reinterpret_cast<uint8x16_t *>(lhs);
  auto * rhs_ptr = reinterpret_cast<uint8x16_t *>(rhs);
  uint64_t count {0};

  for(auto i = 0ULL; i < n_vector_lengths; ++i) {
    count += vaddvq_u8(vcntq_u8(veorq_u8(*lhs_ptr, *rhs_ptr)));
    ++lhs_ptr;
    ++rhs_ptr;
  }

  return count;
}

#elif defined __PPC__

auto compareqgramvectors(unsigned char * lhs, unsigned char * rhs) -> uint64_t
{
  static constexpr auto n_vector_lengths = qgramvectorbytes / sizeof(vector unsigned char);  // 8
  auto * lhs_ptr = reinterpret_cast<vector unsigned char *>(lhs);
  auto * rhs_ptr = reinterpret_cast<vector unsigned char *>(rhs);
  vector unsigned long long count_vector = { 0, 0 };

  for(auto i = 0ULL; i < n_vector_lengths; ++i) {
    count_vector += vec_vpopcnt(reinterpret_cast<vector unsigned long long>(vec_xor(*lhs_ptr, *rhs_ptr)));
    ++lhs_ptr;
    ++rhs_ptr;
  }

  return count_vector[0] + count_vector[1];
}

#elif defined __x86_64__
#ifdef __SSE2__

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


auto compareqgramvectors_128(unsigned char * lhs, unsigned char * rhs) -> uint64_t
{
  /* Count number of different bits */
  /* Uses SSE2 but not POPCNT instruction */
  assert(qgramvectorbytes % 16 == 0); // input MUST be 16-byte aligned

  static constexpr auto n_vector_lengths = qgramvectorbytes / sizeof(__m128i);  // 8
  auto * lhs_ptr = reinterpret_cast<__m128i *>(lhs);
  auto * rhs_ptr = reinterpret_cast<__m128i *>(rhs);
  uint64_t count {0};

  for(auto i = 0ULL; i < n_vector_lengths; ++i) {
    count += popcount_128(_mm_xor_si128(*lhs_ptr, *rhs_ptr));
    ++lhs_ptr;
    ++rhs_ptr;
  }

  return count;
}


auto compareqgramvectors(unsigned char * lhs, unsigned char * rhs) -> uint64_t
{
  if (popcnt_present != 0) {
    return compareqgramvectors_popcnt(lhs, rhs);
  }
  return compareqgramvectors_128(lhs, rhs);
}

#endif
#else

#error Unknown architecture

#endif


inline auto db_getqgramvector(uint64_t seqno) -> unsigned char *
{
  return reinterpret_cast<unsigned char*>(qgrams + seqno);
}


inline auto qgram_diff(uint64_t seqno_a, uint64_t seqno_b) -> uint64_t
{
  const uint64_t diffqgrams = compareqgramvectors(db_getqgramvector(seqno_a),
                                                  db_getqgramvector(seqno_b));
  return (diffqgrams + 2ULL * qgramlength - 1) / (2ULL * qgramlength);  // mindiff
}


auto qgram_worker(int64_t const nth_thread) -> void
{
  auto const & tip = *std::next(thread_info_ptr, nth_thread);

  const auto seed = tip.seed;
  const auto listlen = tip.listlen;
  auto * amplist = tip.amplist;
  auto * difflist = tip.difflist;

  for(auto i = 0ULL; i < listlen; i++) {
    difflist[i] = qgram_diff(seed, amplist[i]);
  }
}


auto qgram_diff_init(std::vector<struct thread_info_s>& thread_info_v) -> void
{
  /* allocate memory for thread info */
  thread_info_v.resize(static_cast<uint64_t>(opt_threads));
  thread_info_ptr = thread_info_v.data();

  qgram_threads
    = new ThreadRunner(static_cast<int>(opt_threads), qgram_worker);
}


auto qgram_diff_done() -> void
{
  delete qgram_threads;
  qgram_threads = nullptr;
  thread_info_ptr = nullptr;
}


auto qgram_diff_fast(uint64_t seed,
                     uint64_t listlen,
                     uint64_t * amplist,
                     uint64_t * difflist,
                     std::vector<struct thread_info_s>& thread_info_v) -> void
{
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
  if (listlen <= uint8_max)
    {
      auto & tip = thread_info_v[0];
      tip.seed = seed;
      tip.listlen = listlen;
      tip.amplist = amplist;
      tip.difflist = difflist;
      qgram_worker(0);
    }
  else
    {
      auto * next_amplist = amplist;
      auto * next_difflist = difflist;
      auto listrest = listlen;
      auto thrrest = static_cast<uint64_t>(opt_threads);

      /* distribute work */
      for(auto & tip: thread_info_v) {
          const uint64_t chunk = (listrest + thrrest - 1) / thrrest;

          tip.seed = seed;
          tip.listlen = chunk;
          tip.amplist = next_amplist;
          tip.difflist = next_difflist;

          next_amplist += chunk;
          next_difflist += chunk;
          listrest -= chunk;
          --thrrest;
        }

      qgram_threads->run();
    }
}
