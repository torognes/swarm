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
#include "utils/qgram_array.h"
#include "utils/threads.h"

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
#endif

#include "popcnt.h"  // refactoring: fence with #ifdef __POPCNT__?

#endif

#include "utils/nt_codec.h"
#include <cassert>
#include <cstdint>  // int64_t, uint64_t
#include <cstring>  // memset


qgramvector_t * qgrams {nullptr};
static ThreadRunner * qgram_threads = nullptr;

static struct thread_info_s
{
  uint64_t seed;
  uint64_t listlen;
  uint64_t * amplist;
  uint64_t * difflist;
} * ti;


void findqgrams(unsigned char * seq, uint64_t seqlen,
                unsigned char * qgramvector)
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

void qgram_work_diff(thread_info_s * tip);
void qgram_worker(int64_t t);
auto compareqgramvectors(unsigned char * qgram_a, unsigned char * qgram_b) -> uint64_t;

#ifdef __aarch64__

uint64_t compareqgramvectors(unsigned char * qgram_a, unsigned char * qgram_b)
{
  uint8x16_t * ap = reinterpret_cast<uint8x16_t *>(qgram_a);
  uint8x16_t * bp = reinterpret_cast<uint8x16_t *>(qgram_b);
  uint64_t count {0};

  while (reinterpret_cast<unsigned char*>(ap) < qgram_a + qgramvectorbytes) {
    count += vaddvq_u8(vcntq_u8(veorq_u8(*ap++, *bp++)));
  }

  return count;
}

#elif defined __PPC__

uint64_t compareqgramvectors(unsigned char * qgram_a, unsigned char * qgram_b)
{
  vector unsigned char * ap = (vector unsigned char *) qgram_a;
  vector unsigned char * bp = (vector unsigned char *) qgram_b;
  vector unsigned long long count_vector = { 0, 0 };

  while ((unsigned char *)ap < qgram_a + qgramvectorbytes) {
    count_vector += vec_vpopcnt((vector unsigned long long)(vec_xor(*ap++, *bp++)));
  }

  return count_vector[0] + count_vector[1];
}

#elif defined __x86_64__

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


auto compareqgramvectors_128(unsigned char * qgram_a, unsigned char * qgram_b) -> uint64_t
{
  /* Count number of different bits */
  /* Uses SSE2 but not POPCNT instruction */
  /* input MUST be 16-byte aligned */

  auto * ap = reinterpret_cast<__m128i *>(qgram_a);
  auto * bp = reinterpret_cast<__m128i *>(qgram_b);
  uint64_t count {0};

  while (reinterpret_cast<unsigned char*>(ap) < qgram_a + qgramvectorbytes) {
    count += popcount_128(_mm_xor_si128(*ap++, *bp++));
  }

  return count;
}


auto compareqgramvectors(unsigned char * qgram_a, unsigned char * qgram_b) -> uint64_t
{
  if (popcnt_present != 0) {
    return compareqgramvectors_popcnt(qgram_a, qgram_b);
  }
  return compareqgramvectors_128(qgram_a, qgram_b);
}

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

auto qgram_worker(int64_t t) -> void
{
  struct thread_info_s * tip = ti + t;

  const uint64_t seed = tip->seed;
  const uint64_t listlen = tip->listlen;
  uint64_t * amplist = tip->amplist;
  uint64_t * difflist = tip->difflist;

  for(auto i = 0ULL; i < listlen; i++) {
    difflist[i] = qgram_diff(seed, amplist[i]);
  }
}

void qgram_diff_init()
{
  /* allocate memory for thread info */
  ti = new struct thread_info_s[static_cast<uint64_t>(opt_threads)];

  qgram_threads
    = new ThreadRunner(static_cast<int> (opt_threads), qgram_worker);
}

void qgram_diff_done()
{
  delete qgram_threads;
  qgram_threads = nullptr;
  delete [] ti;
  ti = nullptr;
}

void qgram_diff_fast(uint64_t seed,
                     uint64_t listlen,
                     uint64_t * amplist,
                     uint64_t * difflist)
{
  if (listlen <= UINT8_MAX)
    {
      ti->seed = seed;
      ti->amplist = amplist;
      ti->difflist = difflist;
      ti->listlen = listlen;
      qgram_worker(0);
    }
  else
    {
      auto thr = static_cast<uint64_t>(opt_threads);

      uint64_t * next_amplist = amplist;
      uint64_t * next_difflist = difflist;
      uint64_t listrest = listlen;
      uint64_t thrrest = thr;

      /* distribute work */
      for(auto t = 0ULL; t < thr; t++)
        {
          thread_info_s * tip = ti + t;
          const uint64_t chunk = (listrest + thrrest - 1) / thrrest;

          tip->seed = seed;
          tip->amplist = next_amplist;
          tip->difflist = next_difflist;
          tip->listlen = chunk;

          next_amplist += chunk;
          next_difflist += chunk;
          listrest -= chunk;
          --thrrest;
        }

      qgram_threads->run();
    }
}
