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

pthread_t pthread_id[MAX_THREADS];

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

    nw(query.seq, query.seq + query.len,
       dseq, dseq + dlen,
       score_matrix_63,
       penalty_gapopen, penalty_gapextend,
       master_scores + sdp->target_index + i,
       master_diffs + sdp->target_index + i,
       master_alignlengths + sdp->target_index + i,
       & nwalignment,
       query.qno, seqno);

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
 
void * worker_8(void * vp)
{
  long t = (long) vp;
  search_init(sd+t);
  while(search_getwork(& sd[t].target_count, & sd[t].target_index))
    search_chunk(sd+t, 8);
  return 0;
}

void * worker_16(void * vp)
{
  long t = (long) vp;
  search_init(sd+t);
  while(search_getwork(& sd[t].target_count, & sd[t].target_index))
    search_chunk(sd+t, 16);
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

  void * (*worker_func)(void *);

  query.qno = query_no;
  db_getsequenceandlength(query_no, &query.seq, &query.len);
  
  master_next = 0;
  master_length = listlength;
  master_targets = targets;
  master_scores = scores;
  master_diffs = diffs;
  master_alignlengths = alignlengths;

  long thr = threads;
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

  if (bits == 16)
    worker_func = worker_16;
  else
    worker_func = worker_8;
  
  if (thr == 1)
  {
    worker_func((void*)0);
  }
  else
  {
    long t;
    void * status;
    
    for(t=0; t<thr; t++)
    {
      
      if (pthread_create(pthread_id + t, 0, worker_func, (void *)t))
	fatal("Cannot create thread.");
    }
    
    for(t=0; t<thr; t++) {
      if (pthread_join(pthread_id[t], &status))
	fatal("Cannot join thread.");
    }
  }
}


void search_begin()
{
  longestdbsequence = db_getlongestsequence();
  
  sd = (struct search_data *) xmalloc(sizeof(search_data) * threads);

  for(int t=0; t<threads; t++)
    search_alloc(sd+t);
}

void search_end()
{
  for(int t=0; t<threads; t++)
    search_free(sd+t);
  free(sd);
}
