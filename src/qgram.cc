/*
    SWARM

    Copyright (C) 2012-2014 Torbjorn Rognes and Frederic Mahe

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

static pthread_attr_t attr;

static struct thread_info_s
{
  /* generic thread info */
  pthread_t pthread;
  pthread_mutex_t workmutex;
  pthread_cond_t workcond;
  int work;

  /* specialized thread info */
  unsigned long seed;
  unsigned long listlen;
  unsigned long * amplist;
  unsigned long * difflist;
} * ti;

void printqgrams(unsigned char * qgramvector)
{
  /* print qgramvector */
  fprintf(logfile, "qgram vector:\n");
  for(int i = 0; i < QGRAMVECTORBYTES; i++)
  {
    fprintf(logfile, "%02x", qgramvector[i]);
    if ((i % 32) == 31)
      fprintf(logfile, "\n");
  }
}

void findqgrams(unsigned char * seq, unsigned long seqlen, 
                unsigned char * qgramvector)
{
  /* set qgram bit vector by xoring occurrences of qgrams in sequence */

  memset(qgramvector, 0, QGRAMVECTORBYTES);
  
  unsigned long qgram = 0;
  unsigned long i = 0;

  while((i < QGRAMLENGTH-1) && (i<seqlen))
  {
    qgram = (qgram << 2) | (seq[i] - 1);
    i++;
  }

  while(i < seqlen)
  {
    qgram = (qgram << 2) | (seq[i] - 1);
    qgramvector[(qgram >> 3) & (QGRAMVECTORBYTES-1)] ^= (1 << (qgram & 7)); 
    i++;
  }
}

/* 
   Unable to get the Mac gcc compiler v 4.2.1 produce the real 
   popcnt instruction. Therefore resorting to assembly code.
*/

#define popcnt_asm(x,y)                                         \
  __asm__ __volatile__ ("popcnt %1,%0" : "=r"(y) : "r"(x));

inline unsigned long popcount(unsigned long x)
{
  unsigned long y;
  popcnt_asm(x,y);
  return y;
}

unsigned long popcount_128(__m128i x)
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

  unsigned long o = (unsigned long) _mm_movepi64_pi64(n);

  return o;
}

unsigned long compareqgramvectors_128(unsigned char * a, unsigned char * b)
{
  /* Count number of different bits */
  /* Uses SSE2 but not POPCNT instruction */
  /* input MUST be 16-byte aligned */

  __m128i * ap = (__m128i *) a;
  __m128i * bp = (__m128i *) b;
  unsigned long count = 0;

  while ((unsigned char*)ap < a + QGRAMVECTORBYTES)
    count += popcount_128(_mm_xor_si128(*ap++, *bp++));
  
  return count;
}

unsigned long compareqgramvectors_64(unsigned char * a, unsigned char * b)
{
  /* Count number of different bits */
  /* Uses the POPCNT instruction, requires CPU with this feature */

  unsigned long *ap = (unsigned long*)a;
  unsigned long *bp = (unsigned long*)b;
  unsigned long count = 0;

  while ((unsigned char*) ap < a + QGRAMVECTORBYTES)
    count += popcount(*ap++ ^ *bp++);
  
  return count;
}

unsigned long compareqgramvectors(unsigned char * a, unsigned char * b)
{
  if (popcnt_present)
    return compareqgramvectors_64(a,b);
  else
    return compareqgramvectors_128(a,b);
}


inline unsigned long qgram_diff(unsigned long a, unsigned long b)
{
  unsigned long diffqgrams = compareqgramvectors(db_getqgramvector(a),
                                                 db_getqgramvector(b));
  unsigned long mindiff = (diffqgrams + 2*QGRAMLENGTH - 1)/(2*QGRAMLENGTH);
  return mindiff;
}

void qgram_work_diff(thread_info_s * tip)
{
  unsigned long seed = tip->seed;
  unsigned long listlen = tip->listlen;
  unsigned long * amplist = tip->amplist;
  unsigned long * difflist = tip->difflist;

  for(unsigned long i=0; i<listlen; i++)
    difflist[i] = qgram_diff(seed, amplist[i]);
}

void * qgram_worker(void * vp)
{
  long t = (long) vp;
  struct thread_info_s * tip = ti + t;

  pthread_mutex_lock(&tip->workmutex);

  /* loop until signalled to quit */
  while (tip->work >= 0)
    {
      /* wait for work available */
      pthread_cond_wait(&tip->workcond, &tip->workmutex);
      if (tip->work > 0)
        {
          qgram_work_diff(tip);
          tip->work = 0;
          pthread_cond_signal(&tip->workcond);
        }
    }
  pthread_mutex_unlock(&tip->workmutex);
  return 0;
}

void qgram_diff_init()
{
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  
  /* allocate memory for thread info */
  ti = (struct thread_info_s *) xmalloc(threads * 
                                        sizeof(struct thread_info_s));
  
  /* init and create worker threads */
  for(unsigned long t=0; t<threads; t++)
    {
      struct thread_info_s * tip = ti + t;
      tip->work = 0;
      pthread_mutex_init(&tip->workmutex, NULL);
      pthread_cond_init(&tip->workcond, NULL);
      if (pthread_create(&tip->pthread, &attr, qgram_worker, (void*)(long)t))
        fatal("Cannot create thread");
    }
}

void qgram_diff_done()
{
  /* finish and clean up worker threads */
  for(unsigned long t=0; t<threads; t++)
    {
      struct thread_info_s * tip = ti + t;
      
      /* tell worker to quit */
      pthread_mutex_lock(&tip->workmutex);
      tip->work = -1;
      pthread_cond_signal(&tip->workcond);
      pthread_mutex_unlock(&tip->workmutex);

      /* wait for worker to quit */
      if (pthread_join(tip->pthread, NULL))
        fatal("Cannot join thread");

      pthread_cond_destroy(&tip->workcond);
      pthread_mutex_destroy(&tip->workmutex);
    }

  free(ti);
  pthread_attr_destroy(&attr);
}

void qgram_diff_fast(unsigned long seed,
                     unsigned long listlen,
                     unsigned long * amplist,
                     unsigned long * difflist)
{
  long thr = threads;
  
  const unsigned long m = 150;

  if (listlen < m*thr)
    thr = (listlen+m-1)/m;
  
  unsigned long * next_amplist = amplist;
  unsigned long * next_difflist = difflist;
  unsigned long listrest = listlen;
  unsigned long thrrest = thr;
  
  /* distribute work */
  for(long t=0; t<thr; t++)
    {
      thread_info_s * tip = ti + t;
      unsigned long chunk = (listrest + thrrest - 1) / thrrest;

      tip->seed = seed;
      tip->amplist = next_amplist;
      tip->difflist = next_difflist;
      tip->listlen = chunk;

      next_amplist += chunk;
      next_difflist += chunk;
      listrest -= chunk;
      thrrest--;
    }

  if (thr == 1)
    {
      qgram_work_diff(ti);
    }
  else
    {
      /* wake up threads */
      for(long t=0; t<thr; t++)
        {
          struct thread_info_s * tip = ti + t;
          pthread_mutex_lock(&tip->workmutex);
          tip->work = 1;
          pthread_cond_signal(&tip->workcond);
          pthread_mutex_unlock(&tip->workmutex);
        }
      
      /* wait for threads to finish their work */
      for(int t=0; t<thr; t++)
        {
          struct thread_info_s * tip = ti + t;
          pthread_mutex_lock(&tip->workmutex);
          while (tip->work > 0)
            pthread_cond_wait(&tip->workcond, &tip->workmutex);
          pthread_mutex_unlock(&tip->workmutex);
        }
    }
}

