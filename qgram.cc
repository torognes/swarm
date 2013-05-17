/*
    SWARM

    Copyright (C) 2012-2013 Torbjorn Rognes and Frederic Mahe

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
#include "limits.h"

unsigned long qgram_map[QGRAMVECTORBITS];

void printqgrams(unsigned char * qgramvector)
{
  /* print qgramvector */
  fprintf(stderr, "qgram vector:\n");
  for(int i = 0; i < QGRAMVECTORBYTES; i++)
  {
    fprintf(stderr, "%02x", qgramvector[i]);
    if ((i % 32) == 31)
      fprintf(stderr, "\n");
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

#define popcnt_asm(x,y)						\
  __asm__ __volatile__ ("popcnt %1,%0" : "=r"(y) : "r"(x));

inline unsigned long popcount(unsigned long x)
{
  unsigned long y;
  popcnt_asm(x,y);
  return y;
}

unsigned long compareqgramvectors(unsigned char * a, unsigned char * b)
{
  /* count number of different bits */
  unsigned long count = 0;

#if 1

  /* using ordinary 64-bit registers */

  unsigned long *ap = (unsigned long*)a;
  unsigned long *bp = (unsigned long*)b;

#if 1

  /* with ordinary loop */

  for(int i = 0; i < QGRAMVECTORBYTES/8; i++)
    count += popcount(ap[i] ^ bp[i]);

#else

  /* unrolled loop */

  count += popcount(ap[ 0] ^ bp[ 0]);
  count += popcount(ap[ 1] ^ bp[ 1]);
  count += popcount(ap[ 2] ^ bp[ 2]);
  count += popcount(ap[ 3] ^ bp[ 3]);
  count += popcount(ap[ 4] ^ bp[ 4]);
  count += popcount(ap[ 5] ^ bp[ 5]);
  count += popcount(ap[ 6] ^ bp[ 6]);
  count += popcount(ap[ 7] ^ bp[ 7]);
  count += popcount(ap[ 8] ^ bp[ 8]);
  count += popcount(ap[ 9] ^ bp[ 9]);
  count += popcount(ap[10] ^ bp[10]);
  count += popcount(ap[11] ^ bp[11]);
  count += popcount(ap[12] ^ bp[12]);
  count += popcount(ap[13] ^ bp[13]);
  count += popcount(ap[14] ^ bp[14]);
  count += popcount(ap[15] ^ bp[15]);

#endif

#else

  /* SSE2 version using the 128-bit XMM registers */

  __m128i * ap = (__m128i*) a;
  __m128i * bp = (__m128i*) b;

  for(int i = 0; i < 8; i++)
    {
      __m128i x = _mm_load_si128(ap+i);
      x = _mm_xor_si128(x, bp[i]);
      unsigned long u = _mm_extract_epi64(x, 0);
      unsigned long v = _mm_extract_epi64(x, 1);
      count += popcount(u);
      count += popcount(v);
    }

#endif

  return count;
}

inline unsigned long qgram_diff(unsigned long a, unsigned long b)
{
  unsigned long diffqgrams = compareqgramvectors(db_getqgramvector(a),
						 db_getqgramvector(b));
  unsigned long mindiff = (diffqgrams + 2*QGRAMLENGTH - 1)/(2*QGRAMLENGTH);
  return mindiff;
}

struct thread_info_struct
{
  pthread_t id;
  unsigned long seed;
  unsigned long listlen;
  unsigned long * amplist;
  unsigned long * difflist;
} thread_info[MAX_THREADS];

void * qgram_diff_worker(void * vp)
{
  long t = (long) vp;

  unsigned long seed = thread_info[t].seed;
  unsigned long listlen = thread_info[t].listlen;
  unsigned long * amplist = thread_info[t].amplist;
  unsigned long * difflist = thread_info[t].difflist;

  for(unsigned long i=0; i<listlen; i++)
    difflist[i] = qgram_diff(seed, amplist[i]);
  return 0;
}

void qgram_diff_parallel(unsigned long seed,
			 unsigned long listlen,
			 unsigned long * amplist,
			 unsigned long * difflist)
{
  long thr = threads;
  
  const unsigned long m = 3000;

  if (listlen < m*thr)
    thr = (listlen+m-1)/m;
  
  if (thr == 1)
    {
      for(unsigned long i=0; i<listlen; i++)
	difflist[i] = qgram_diff(seed, amplist[i]);
    }
  else
    {
      long t;
      void * status;
      
      unsigned long * next_amplist = amplist;
      unsigned long * next_difflist = difflist;
      unsigned long listrest = listlen;
      unsigned long thrrest = thr;
      unsigned long chunk;
      
      for(t=0; t<thr; t++)
	{
	  thread_info[t].seed = seed;
	  thread_info[t].amplist = next_amplist;
	  thread_info[t].difflist = next_difflist;
	  chunk = (listrest + thrrest - 1) / thrrest;
	  thread_info[t].listlen = chunk;
	  next_amplist += chunk;
	  next_difflist += chunk;
	  listrest -= chunk;
	  thrrest--;
	}

      for(t=0; t<thr; t++)
	if (pthread_create(&thread_info[t].id, 0, qgram_diff_worker, (void *)t))
	  fatal("Cannot create thread.");
      
      for(t=0; t<thr; t++) {
	if (pthread_join(thread_info[t].id, &status))
	  fatal("Cannot join thread.");
      }
    }
}
