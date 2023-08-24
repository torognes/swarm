/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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
#include "threads.h"
#include "util.h"

#ifdef __x86_64__

#include "popcnt.h"

#endif



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

  memset(qgramvector, 0, qgramvectorbytes);

  uint64_t qgram {0};
  unsigned int i {0};

  while((i < qgramlength-1) && (i<seqlen))
  {
    qgram = (qgram << 2U) | nt_extract(reinterpret_cast<char *>(seq), i);
    i++;
  }

  while(i < seqlen)
  {
    qgram = (qgram << 2U) | nt_extract(reinterpret_cast<char *>(seq), i);
    qgramvector[(qgram >> 3U) & (qgramvectorbytes - 1)] ^= (1U << (qgram & max_range));
    i++;
  }
}

void qgram_work_diff(thread_info_s * tip);
void qgram_worker(int64_t t);
auto compareqgramvectors(unsigned char * a, unsigned char * b) -> uint64_t;

#ifdef __aarch64__

uint64_t compareqgramvectors(unsigned char * a, unsigned char * b)
{
  uint8x16_t * ap = (uint8x16_t *) a;
  uint8x16_t * bp = (uint8x16_t *) b;
  uint64_t count {0};

  while ((unsigned char*)ap < a + qgramvectorbytes) {
    count += vaddvq_u8(vcntq_u8(veorq_u8(*ap++, *bp++)));
  }

  return count;
}

#elif defined __PPC__

uint64_t compareqgramvectors(unsigned char * a, unsigned char * b)
{
  vector unsigned char * ap = (vector unsigned char *) a;
  vector unsigned char * bp = (vector unsigned char *) b;
  vector unsigned long long count_vector = { 0, 0 };

  while ((unsigned char *)ap < a + qgramvectorbytes) {
    count_vector += vec_vpopcnt((vector unsigned long long)(vec_xor(*ap++, *bp++)));
  }

  return count_vector[0] + count_vector[1];
}

#elif defined __x86_64__

auto popcount_128(__m128i x) -> uint64_t;
auto compareqgramvectors_128(unsigned char * a, unsigned char * b) -> uint64_t;

auto popcount_128(__m128i x) -> uint64_t
{
  static constexpr unsigned char m1 {0x55};  // '0101 0101'
  static constexpr unsigned char m2 {0x33};  // '0011 0011'
  static constexpr unsigned char m4 {0x0f};  // '0000 1111'

  const __m128i mask1 = _mm_set_epi8(m1, m1, m1, m1, m1, m1, m1, m1,
                                     m1, m1, m1, m1, m1, m1, m1, m1);

  const __m128i mask2 = _mm_set_epi8(m2, m2, m2, m2, m2, m2, m2, m2,
                                     m2, m2, m2, m2, m2, m2, m2, m2);

  const __m128i mask4 = _mm_set_epi8(m4, m4, m4, m4, m4, m4, m4, m4,
                                     m4, m4, m4, m4, m4, m4, m4, m4);

  const __m128i zero = _mm_setzero_si128();

  /* add together 2 bits: 0+1, 2+3, 3+4, ... 126+127 */

  const __m128i a = _mm_srli_epi64(x, 1);
  const __m128i b = _mm_and_si128(x, mask1);
  const __m128i c = _mm_and_si128(a, mask1);
  const __m128i d = _mm_add_epi64(b, c);

  /* add together 4 bits: (0+1)+(2+3), ... (124+125)+(126+127) */

  const __m128i e = _mm_srli_epi64(d, 2);
  const __m128i f = _mm_and_si128(d, mask2);
  const __m128i g = _mm_and_si128(e, mask2);
  const __m128i h = _mm_add_epi64(f, g);

  /* add together 8 bits: (0..3)+(4..7), ... (120..123)+(124..127) */

  const __m128i i = _mm_srli_epi64(h, 4);
  const __m128i j = _mm_add_epi64(h, i);
  const __m128i k = _mm_and_si128(j, mask4);

  /* add together 8 bytes: (0..63) and (64..127) */

  const __m128i l = _mm_sad_epu8(k, zero);

  /* add together 64-bit values into final 128 bit value */

  const __m128i m = _mm_srli_si128(l, 8);
  const __m128i n = _mm_add_epi64(m, l);

  /* return low 64 bits: return value is always in range 0 to 128 */

  return reinterpret_cast<uint64_t>(_mm_movepi64_pi64(n));
}

auto compareqgramvectors_128(unsigned char * a, unsigned char * b) -> uint64_t
{
  /* Count number of different bits */
  /* Uses SSE2 but not POPCNT instruction */
  /* input MUST be 16-byte aligned */

  auto * ap = reinterpret_cast<__m128i *>(a);
  auto * bp = reinterpret_cast<__m128i *>(b);
  uint64_t count {0};

  while (reinterpret_cast<unsigned char*>(ap) < a + qgramvectorbytes) {
    count += popcount_128(_mm_xor_si128(*ap++, *bp++));
  }

  return count;
}

auto compareqgramvectors(unsigned char * a, unsigned char * b) -> uint64_t
{
  if (popcnt_present != 0) {
    return compareqgramvectors_popcnt(a, b);
  }
  return compareqgramvectors_128(a, b);
}

#else

#error Unknown architecture

#endif

inline auto db_getqgramvector(uint64_t seqno) -> unsigned char *
{
  return reinterpret_cast<unsigned char*>(qgrams + seqno);
}

inline auto qgram_diff(uint64_t a, uint64_t b) -> uint64_t
{
  const uint64_t diffqgrams = compareqgramvectors(db_getqgramvector(a),
                                                  db_getqgramvector(b));
  return (diffqgrams + 2 * qgramlength - 1)/(2 * qgramlength);  // mindiff
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
          thrrest--;
        }

      qgram_threads->run();
    }
}
