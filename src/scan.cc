/*
    SWARM

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

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


pthread_mutex_t workmutex = PTHREAD_MUTEX_INITIALIZER;

queryinfo_t query;

struct search_data
{
  BYTE ** qtable;
  WORD ** qtable_w;

  BYTE * dprofile;
  WORD * dprofile_w;

  BYTE * hearray;

  unsigned long * dir_array;

  unsigned long target_count;
  unsigned long target_index;
};

struct search_data * sd;

unsigned long master_next;
unsigned long master_length;

unsigned long remainingchunks;

unsigned long * master_targets;
unsigned long * master_scores;
unsigned long * master_diffs;
unsigned long * master_alignlengths;
int master_bits;

unsigned long longestdbsequence;
unsigned long dirbufferbytes;

void search_alloc(struct search_data * sdp)
{
  dirbufferbytes = 8 * longestdbsequence * ((longestdbsequence+3)/4) * 4;
  sdp->qtable = (BYTE**) xmalloc(longestdbsequence * sizeof(BYTE*));
  sdp->qtable_w = (WORD**) xmalloc(longestdbsequence * sizeof(WORD*));
  sdp->dprofile = (BYTE*) xmalloc(4*16*32);
  sdp->dprofile_w = (WORD*) xmalloc(4*2*8*32);
  sdp->hearray = (BYTE*) xmalloc(longestdbsequence * 32);
  sdp->dir_array = (unsigned long *) xmalloc(dirbufferbytes);

  memset(sdp->hearray, 0, longestdbsequence*32);
  memset(sdp->dir_array, 0, dirbufferbytes);
}

void search_free(struct search_data * sdp)
{
  free(sdp->qtable);
  free(sdp->qtable_w);
  free(sdp->dprofile);
  free(sdp->dprofile_w);
  free(sdp->hearray);
  free(sdp->dir_array);
}

void search_init(struct search_data * sdp)
{
  for (long i = 0; i < query.len; i++ )
  {
    sdp->qtable[i] = sdp->dprofile + 64*query.seq[i];
    sdp->qtable_w[i] = sdp->dprofile_w + 32*query.seq[i];
  }
}

void search_chunk(struct search_data * sdp, long bits)
{
  if (sdp->target_count == 0)
    return;

#if 0

  for(unsigned long i=0; i<sdp->target_count; i++)
    {
    char * dseq;
    long dlen;
    char * nwalignment;

    unsigned long seqno = master_targets[sdp->target_index + i];
    db_getsequenceandlength(seqno, & dseq, & dlen);

    nw(dseq, dseq + dlen,
       query.seq, query.seq + query.len,
       score_matrix_63,
       penalty_gapopen, penalty_gapextend,
       master_scores + sdp->target_index + i,
       master_diffs + sdp->target_index + i,
       master_alignlengths + sdp->target_index + i,
       & nwalignment,
       (unsigned char *) sdp->dir_array,
       (unsigned long int *) sdp->hearray,
       query.qno, seqno);

#if 0
    printf("\nAlignment: %s\n", nwalignment);
#endif

    free(nwalignment);
  }

  return;

#endif

  if (bits == 16)
    search16(sdp->qtable_w,
             penalty_gapopen,
             penalty_gapextend,
             (WORD*) score_matrix_16,
             sdp->dprofile_w,
             (WORD*) sdp->hearray,
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
            penalty_gapopen,
            penalty_gapextend,
            (BYTE*) score_matrix_8,
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

int search_getwork(unsigned long * countref, unsigned long * firstref)
{
  // * countref = how many sequences to search
  // * firstref = index into master_targets/scores/diffs where thread should start

  unsigned long status = 0;

  pthread_mutex_lock(&workmutex);

  if (master_next < master_length)
    {
      unsigned long chunksize =
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

void master_dump()
{
  printf("master_dump\n");
  printf("   i    t    s    d\n");
  for(unsigned long i=0; i< 1403; i++)
    {
      printf("%4lu %4lu %4lu %4lu\n", i, master_targets[i],
             master_scores[i], master_diffs[i]);
    }
}

void search_worker_core(int t)
{
  search_init(sd+t);
  while(search_getwork(& sd[t].target_count, & sd[t].target_index))
    search_chunk(sd+t, master_bits);
}

void * search_worker(void * vp)
{
  long t = (long) vp;
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
  return 0;
}

void search_do(unsigned long query_no,
               unsigned long listlength,
               unsigned long * targets,
               unsigned long * scores,
               unsigned long * diffs,
               unsigned long * alignlengths,
               long bits)
{
  query.qno = query_no;
  db_getsequenceandlength(query_no, &query.seq, &query.len);

  master_next = 0;
  master_length = listlength;
  master_targets = targets;
  master_scores = scores;
  master_diffs = diffs;
  master_alignlengths = alignlengths;
  master_bits = bits;

  unsigned long thr = opt_threads;

  if (bits == 8)
    {
      if (master_length <= (unsigned long)(15 * thr) )
        thr = (master_length + 15) / 16;
    }
  else
    {
      if (master_length <= (unsigned long)(7 * thr) )
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
      for(unsigned long t=0; t<thr; t++)
        {
          struct thread_info_s * tip = ti + t;
          pthread_mutex_lock(&tip->workmutex);
          tip->work = 1;
          pthread_cond_signal(&tip->workcond);
          pthread_mutex_unlock(&tip->workmutex);
        }

      /* wait for threads to finish their work */
      for(unsigned long t=0; t<thr; t++)
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

  sd = (struct search_data *) xmalloc(sizeof(search_data) * opt_threads);

  for(long t=0; t<opt_threads; t++)
    search_alloc(sd+t);

  /* start threads */

  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* allocate memory for thread info */
  ti = (struct thread_info_s *) xmalloc(opt_threads *
                                        sizeof(struct thread_info_s));

  /* init and create worker threads */
  for(long t=0; t<opt_threads; t++)
    {
      struct thread_info_s * tip = ti + t;
      tip->work = 0;
      pthread_mutex_init(&tip->workmutex, NULL);
      pthread_cond_init(&tip->workcond, NULL);
      if (pthread_create(&tip->pthread, &attr, search_worker, (void*)(long)t))
        fatal("Cannot create thread");
    }

}

void search_end()
{
  /* finish and clean up worker threads */

  for(long t=0; t<opt_threads; t++)
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

  for(long t=0; t<opt_threads; t++)
    search_free(sd+t);
  free(sd);
}
