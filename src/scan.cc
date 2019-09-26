/*
    SWARM

    Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe

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
  uint64_t seed;
  uint64_t listlen;
  uint64_t * amplist;
  uint64_t * difflist;
} * ti;


static pthread_mutex_t workmutex = PTHREAD_MUTEX_INITIALIZER;


struct search_data
{
  BYTE ** qtable;
  WORD ** qtable_w;

  BYTE * dprofile;
  WORD * dprofile_w;

  BYTE * hearray;

  uint64_t * dir_array;

  uint64_t target_count;
  uint64_t target_index;
};

static struct search_data * sd;
static uint64_t master_next;
static uint64_t master_length;
static uint64_t remainingchunks;
static uint64_t * master_targets;
static uint64_t * master_scores;
static uint64_t * master_diffs;
static uint64_t * master_alignlengths;
static int master_bits;
static uint64_t dirbufferbytes;

queryinfo_t query;
uint64_t longestdbsequence;

void search_alloc(struct search_data * sdp)
{
  dirbufferbytes = 8 * longestdbsequence * ((longestdbsequence+3)/4) * 4;
  sdp->qtable = static_cast<BYTE**>
    (xmalloc(longestdbsequence * sizeof(BYTE*)));
  sdp->qtable_w = static_cast<WORD**>
    (xmalloc(longestdbsequence * sizeof(WORD*)));
  sdp->dprofile = static_cast<BYTE*>
    (xmalloc(4*16*32));
  sdp->dprofile_w = static_cast<WORD*>
    (xmalloc(4*2*8*32));
  sdp->hearray = static_cast<BYTE*>
    (xmalloc(longestdbsequence * 32));
  sdp->dir_array = static_cast<uint64_t *>
    (xmalloc(dirbufferbytes));

  memset(sdp->hearray, 0, longestdbsequence*32);
  memset(sdp->dir_array, 0, dirbufferbytes);
}

void search_free(struct search_data * sdp)
{
  xfree(sdp->qtable);
  xfree(sdp->qtable_w);
  xfree(sdp->dprofile);
  xfree(sdp->dprofile_w);
  xfree(sdp->hearray);
  xfree(sdp->dir_array);
}

void search_init(struct search_data * sdp)
{
  for (unsigned int i = 0; i < query.len; i++ )
  {
    sdp->qtable[i] = sdp->dprofile + 64 * (nt_extract(query.seq, i) + 1);
    sdp->qtable_w[i] = sdp->dprofile_w + 32 * (nt_extract(query.seq, i) + 1);
  }
}

void search_chunk(struct search_data * sdp, int64_t bits)
{
  if (sdp->target_count == 0)
    return;

#if 0

  for(uint64_t i=0; i<sdp->target_count; i++)
    {
    char * dseq;
    int64_t dlen;
    char * nwalignment;

    uint64_t seqno = master_targets[sdp->target_index + i];
    db_getsequenceandlength(seqno, & dseq, & dlen);

    nw(dseq, dlen,
       query.seq, query.len,
       score_matrix_63,
       penalty_gapopen, penalty_gapextend,
       master_scores + sdp->target_index + i,
       master_diffs + sdp->target_index + i,
       master_alignlengths + sdp->target_index + i,
       & nwalignment,
       (unsigned char *) sdp->dir_array,
       (uint64_t int *) sdp->hearray,
       query.qno, seqno);

#if 0
    printf("\nAlignment: %s\n", nwalignment);
#endif

    xfree(nwalignment);
  }

  return;

#endif

#ifdef __aarch64__
  /* always use 16-bit version on aarch64 because it is faster */
 (void) bits;
  if (1)
#else
  if (bits == 16)
#endif
    search16(sdp->qtable_w,
             static_cast<WORD>(penalty_gapopen),
             static_cast<WORD>(penalty_gapextend),
             static_cast<WORD*>(score_matrix_16),
             sdp->dprofile_w,
             reinterpret_cast<WORD*>(sdp->hearray),
             sdp->target_count,
             master_targets + sdp->target_index,
             master_scores + sdp->target_index,
             master_diffs + sdp->target_index,
             master_alignlengths + sdp->target_index,
             query.len,
             dirbufferbytes/8,
             sdp->dir_array);
  else
    search8(sdp->qtable,
            static_cast<BYTE>(penalty_gapopen),
            static_cast<BYTE>(penalty_gapextend),
            static_cast<BYTE*>(score_matrix_8),
            sdp->dprofile,
            sdp->hearray,
            sdp->target_count,
            master_targets + sdp->target_index,
            master_scores + sdp->target_index,
            master_diffs + sdp->target_index,
            master_alignlengths + sdp->target_index,
            query.len,
            dirbufferbytes/8,
            sdp->dir_array);
}

int search_getwork(uint64_t * countref, uint64_t * firstref)
{
  // * countref = how many sequences to search
  // * firstref = index into master_targets/scores/diffs where thread should start

  int status = 0;

  pthread_mutex_lock(&workmutex);

  if (master_next < master_length)
    {
      uint64_t chunksize =
        ((master_length - master_next + remainingchunks - 1) / remainingchunks);

      * countref = chunksize;
      * firstref = master_next;

      master_next += chunksize;
      remainingchunks--;
      status = 1;
    }

  pthread_mutex_unlock(&workmutex);

  return status;
}

#if 0

/* never used */

void master_dump()
{
  printf("master_dump\n");
  printf("   i    t    s    d\n");
  for(uint64_t i=0; i< 1403; i++)
    {
      printf("%4" PRIu64 " %4" PRIu64 " %4" PRIu64 " %4" PRIu64 "\n",
             i, master_targets[i], master_scores[i], master_diffs[i]);
    }
}

#endif

void search_worker_core(uint64_t t)
{
  search_init(sd+t);
  while(search_getwork(& sd[t].target_count, & sd[t].target_index))
    search_chunk(sd+t, master_bits);
}

void * search_worker(void * vp)
{
  uint64_t t = reinterpret_cast<uint64_t>(vp);
  struct thread_info_s * tip = ti + t;

  pthread_mutex_lock(&tip->workmutex);

  /* loop until signalled to quit */
  while (tip->work >= 0)
    {
      /* wait for work available */
      while (tip->work == 0)
        pthread_cond_wait(&tip->workcond, &tip->workmutex);
      if (tip->work > 0)
        {
          search_worker_core(t);
          tip->work = 0;
          pthread_cond_signal(&tip->workcond);
        }
    }
  pthread_mutex_unlock(&tip->workmutex);
  return nullptr;
}

void search_do(uint64_t query_no,
               uint64_t listlength,
               uint64_t * targets,
               uint64_t * scores,
               uint64_t * diffs,
               uint64_t * alignlengths,
               int bits)
{
  query.qno = query_no;
  unsigned int query_len = 0;
  db_getsequenceandlength(query_no, &query.seq, &query_len);
  query.len = query_len;

  master_next = 0;
  master_length = listlength;
  master_targets = targets;
  master_scores = scores;
  master_diffs = diffs;
  master_alignlengths = alignlengths;
  master_bits = bits;

  uint64_t thr = opt_threads;

  if (bits == 8)
    {
      if (master_length <= 15 * thr)
        thr = (master_length + 15) / 16;
    }
  else
    {
      if (master_length <= 7 * thr)
        thr = (master_length + 7) / 8;
    }

  remainingchunks = thr;

  if (thr == 1)
    {
      search_worker_core(0);
    }
  else
    {
      /* wake up threads */
      for(uint64_t t=0; t<thr; t++)
        {
          struct thread_info_s * tip = ti + t;
          pthread_mutex_lock(&tip->workmutex);
          tip->work = 1;
          pthread_cond_signal(&tip->workcond);
          pthread_mutex_unlock(&tip->workmutex);
        }

      /* wait for threads to finish their work */
      for(uint64_t t=0; t<thr; t++)
        {
          struct thread_info_s * tip = ti + t;
          pthread_mutex_lock(&tip->workmutex);
          while (tip->work > 0)
            pthread_cond_wait(&tip->workcond, &tip->workmutex);
          pthread_mutex_unlock(&tip->workmutex);
        }
    }
}

void search_begin()
{
  longestdbsequence = db_getlongestsequence();

  sd = static_cast<struct search_data *>
    (xmalloc(sizeof(search_data) * opt_threads));

  for(int64_t t=0; t<opt_threads; t++)
    search_alloc(sd+t);

  /* start threads */

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* allocate memory for thread info */
  ti = static_cast<struct thread_info_s *>
    (xmalloc(opt_threads * sizeof(struct thread_info_s)));

  /* init and create worker threads */
  for(int64_t t=0; t<opt_threads; t++)
    {
      struct thread_info_s * tip = ti + t;
      tip->work = 0;
      pthread_mutex_init(&tip->workmutex, nullptr);
      pthread_cond_init(&tip->workcond, nullptr);
      if (pthread_create(&tip->pthread,
                         &attr,
                         search_worker,
                         reinterpret_cast<void*>(t)))
        fatal("Cannot create thread");
    }
}

void search_end()
{
  /* finish and clean up worker threads */

  for(int64_t t=0; t<opt_threads; t++)
    {
      struct thread_info_s * tip = ti + t;

      /* tell worker to quit */
      pthread_mutex_lock(&tip->workmutex);
      tip->work = -1;
      pthread_cond_signal(&tip->workcond);
      pthread_mutex_unlock(&tip->workmutex);

      /* wait for worker to quit */
      if (pthread_join(tip->pthread, nullptr))
        fatal("Cannot join thread");

      pthread_cond_destroy(&tip->workcond);
      pthread_mutex_destroy(&tip->workmutex);
    }

  xfree(ti);
  pthread_attr_destroy(&attr);

  for(int64_t t=0; t<opt_threads; t++)
    search_free(sd+t);
  xfree(sd);
}
