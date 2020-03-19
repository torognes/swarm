/*
    SWARM

    Copyright (C) 2012-2020 Torbjorn Rognes and Frederic Mahe

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

static ThreadRunner * qgram_threads = nullptr;

static struct thread_info_s
{
  uint64_t seed;
  uint64_t listlen;
  uint64_t * amplist;
  uint64_t * difflist;
} * ti;

#if 0

/* never used */

void printqgrams(unsigned char * qgramvector)
{
  /* print qgramvector */
  fprintf(logfile, "qgram vector:\n");
  for(auto i = 0; i < qgramvectorbytes; i++)
  {
    fprintf(logfile, "%02x", qgramvector[i]);
    if ((i % 32) == 31)
      fprintf(logfile, "\n");
  }
}

#endif

void findqgrams(unsigned char * seq, uint64_t seqlen,
                unsigned char * qgramvector)
{
  /* set qgram bit vector by xoring occurrences of qgrams in sequence */

  memset(qgramvector, 0, qgramvectorbytes);

  uint64_t qgram = 0;
  unsigned int i = 0;

  while((i < qgramlength-1) && (i<seqlen))
  {
    qgram = (qgram << 2) | nt_extract(reinterpret_cast<char *>(seq), i);
    i++;
  }

  while(i < seqlen)
  {
    qgram = (qgram << 2) | nt_extract(reinterpret_cast<char *>(seq), i);
    qgramvector[(qgram >> 3) & (qgramvectorbytes-1)] ^= (1 << (qgram & 7));
    i++;
  }
}

void qgram_work_diff(thread_info_s * tip);
void qgram_worker(int64_t t);
uint64_t compareqgramvectors(unsigned char * a, unsigned char * b);

#ifdef __aarch64__

uint64_t compareqgramvectors(unsigned char * a, unsigned char * b)
{
  uint8x16_t * ap = (uint8x16_t *) a;
  uint8x16_t * bp = (uint8x16_t *) b;
  uint64_t count = 0;

  while ((unsigned char*)ap < a + qgramvectorbytes)
    count += vaddvq_u8(vcntq_u8(veorq_u8(*ap++, *bp++)));

  return count;
}

#elif defined __PPC__

uint64_t compareqgramvectors(unsigned char * a, unsigned char * b)
{
  vector unsigned char * ap = (vector unsigned char *) a;
  vector unsigned char * bp = (vector unsigned char *) b;
  vector unsigned long long count_vector = { 0, 0 };

  while ((unsigned char *)ap < a + qgramvectorbytes)
    count_vector += vec_vpopcnt((vector unsigned long long)(vec_xor(*ap++, *bp++)));

  return count_vector[0] + count_vector[1];
}

#elif defined __x86_64__

/*
   Unable to get the Mac gcc compiler v 4.2.1 produce the real
   popcnt instruction. Therefore resorting to assembly code.
*/

#define popcnt_asm(x,y)                                         \
  __asm__ __volatile__ ("popcnt %1,%0" : "=r"(y) : "r"(x))

inline uint64_t popcount(uint64_t x)
{
  uint64_t y;
  popcnt_asm(x,y);
  return y;
}

uint64_t popcount_128(__m128i x);
uint64_t compareqgramvectors_128(unsigned char * a, unsigned char * b);
uint64_t compareqgramvectors_64(unsigned char * a, unsigned char * b);

uint64_t popcount_128(__m128i x)
{
  __m128i mask1 = _mm_set_epi8(0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55,
                               0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55, 0x55);

  __m128i mask2 = _mm_set_epi8(0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33,
                               0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33, 0x33);

  __m128i mask4 = _mm_set_epi8(0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f,
                               0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f, 0x0f);

  __m128i zero = _mm_setzero_si128();

  /* add together 2 bits: 0+1, 2+3, 3+4, ... 126+127 */

  __m128i a = _mm_srli_epi64(x, 1);
  __m128i b = _mm_and_si128(x, mask1);
  __m128i c = _mm_and_si128(a, mask1);
  __m128i d = _mm_add_epi64(b, c);

  /* add together 4 bits: (0+1)+(2+3), ... (124+125)+(126+127) */

  __m128i e = _mm_srli_epi64(d, 2);
  __m128i f = _mm_and_si128(d, mask2);
  __m128i g = _mm_and_si128(e, mask2);
  __m128i h = _mm_add_epi64(f, g);

  /* add together 8 bits: (0..3)+(4..7), ... (120..123)+(124..127) */

  __m128i i = _mm_srli_epi64(h, 4);
  __m128i j = _mm_add_epi64(h, i);
  __m128i k = _mm_and_si128(j, mask4);

  /* add together 8 bytes: (0..63) and (64..127) */

  __m128i l = _mm_sad_epu8(k, zero);

  /* add together 64-bit values into final 128 bit value */

  __m128i m = _mm_srli_si128(l, 8);
  __m128i n = _mm_add_epi64(m, l);

  /* return low 64 bits: return value is always in range 0 to 128 */

  uint64_t o = reinterpret_cast<uint64_t>(_mm_movepi64_pi64(n));

  return o;
}

uint64_t compareqgramvectors_128(unsigned char * a, unsigned char * b)
{
  /* Count number of different bits */
  /* Uses SSE2 but not POPCNT instruction */
  /* input MUST be 16-byte aligned */

  __m128i * ap = reinterpret_cast<__m128i *>(a);
  __m128i * bp = reinterpret_cast<__m128i *>(b);
  uint64_t count = 0;

  while (reinterpret_cast<unsigned char*>(ap) < a + qgramvectorbytes)
    count += popcount_128(_mm_xor_si128(*ap++, *bp++));

  return count;
}


uint64_t compareqgramvectors_64(unsigned char * a, unsigned char * b)
{
  /* Count number of different bits */
  /* Uses the POPCNT instruction, requires CPU with this feature */

  uint64_t *ap = reinterpret_cast<uint64_t*>(a);
  uint64_t *bp = reinterpret_cast<uint64_t*>(b);
  uint64_t count = 0;

  while (reinterpret_cast<unsigned char*>(ap) < a + qgramvectorbytes)
    count += popcount(*ap++ ^ *bp++);

  return count;
}


uint64_t compareqgramvectors(unsigned char * a, unsigned char * b)
{
  if (popcnt_present)
    return compareqgramvectors_64(a,b);
  else
    return compareqgramvectors_128(a,b);
}

#else

#error Unknown architecture

#endif

inline uint64_t qgram_diff(uint64_t a, uint64_t b)
{
  uint64_t diffqgrams = compareqgramvectors(db_getqgramvector(a),
                                            db_getqgramvector(b));
  uint64_t mindiff = (diffqgrams + 2*qgramlength - 1)/(2*qgramlength);
  return mindiff;
}

void qgram_worker(int64_t t)
{
  struct thread_info_s * tip = ti + t;

  uint64_t seed = tip->seed;
  uint64_t listlen = tip->listlen;
  uint64_t * amplist = tip->amplist;
  uint64_t * difflist = tip->difflist;

  for(auto i = 0ULL; i < listlen; i++)
    difflist[i] = qgram_diff(seed, amplist[i]);
}

void qgram_diff_init()
{
  /* allocate memory for thread info */
  ti = static_cast<struct thread_info_s *>
    (xmalloc(static_cast<uint64_t>(opt_threads) *
             sizeof(struct thread_info_s)));

  qgram_threads
    = new ThreadRunner(static_cast<int> (opt_threads), qgram_worker);
}

void qgram_diff_done()
{
  delete qgram_threads;

  xfree(ti);
}

void qgram_diff_fast(uint64_t seed,
                     uint64_t listlen,
                     uint64_t * amplist,
                     uint64_t * difflist)
{
  if (listlen < 256)
    {
      ti->seed = seed;
      ti->amplist = amplist;
      ti->difflist = difflist;
      ti->listlen = listlen;
      qgram_worker(0);
    }
  else
    {
      uint64_t thr = static_cast<uint64_t>(opt_threads);

      uint64_t * next_amplist = amplist;
      uint64_t * next_difflist = difflist;
      uint64_t listrest = listlen;
      uint64_t thrrest = thr;

      /* distribute work */
      for(auto t = 0ULL; t < thr; t++)
        {
          thread_info_s * tip = ti + t;
          uint64_t chunk = (listrest + thrrest - 1) / thrrest;

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
