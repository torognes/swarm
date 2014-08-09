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

/*
  This version of the swarm algorithm uses Frederic's idea for d=1 to
  enumerate all of the maximum 7L+4 possible variants of a sequence with only
  one difference, where L is the length of the sequence.
*/

#include "swarm.h"

#define SEPCHAR ' '
#define HASH hash_djb2
//#define HASH hash_djb2a
//#define HASH hash_fnv_1a_32
//#define HASH hash_fnv_1a_64

#define HASHFILLFACTOR 4

struct bucket_s
{
  unsigned long hash;
  int seqno;
  int swarmid;
  int generation;
  struct bucket_s * swarms_next;
  struct bucket_s * swarm_next;
  struct bucket_s * all_next;
};

struct bucket_s * all_head;
struct bucket_s * all_tail;
struct bucket_s * swarms_head;
struct bucket_s * swarms_tail;
struct bucket_s * seed;
struct bucket_s * current_swarm_tail;

unsigned long amphashtablesize = 0;
struct bucket_s * amphashtable = 0;

/* overall statistics */
static unsigned long maxgen = 0;
static unsigned long largest = 0;

/* per swarm statistics */
static unsigned long singletons = 0;
static unsigned long abundance_sum = 0;
static unsigned long swarmsize = 0;
static unsigned long swarm_maxgen = 0;

pthread_attr_t attr;
pthread_mutex_t mutex_varmatch;

static struct thread_info_s
{
  pthread_t pthread;
  pthread_mutex_t workmutex;
  pthread_cond_t workcond;
  int work;
  unsigned char * varseq;
  struct bucket_s * seed;
  unsigned long mut_start;
  unsigned long mut_length;
} * ti;

long lastseed = -1;
long seedmatches = 0;

void find_variant_matches(unsigned long thread,
			  unsigned char * seq,
			  unsigned long seqlen,
			  struct bucket_s * seed)
{
  /* compute hash and corresponding hash table index */

  unsigned long hash = HASH(seq, seqlen);
  unsigned long j = hash % amphashtablesize;
  
  /* find matching buckets */

  while (amphashtable[j].seqno >= 0)
    {
      struct bucket_s * bp = amphashtable + j;
      if (bp->hash == hash)
	{
	  unsigned long seqno = bp->seqno;
	  unsigned long ampseqlen = db_getsequencelen(seqno);
	  unsigned char * ampseq = (unsigned char *) db_getsequence(seqno);

	  /* check if not already swarmed */
	  /* make sure sequences are identical even though hashes are */

	  if ((!bp->swarmid) &&
	      (ampseqlen == seqlen) &&
	      (!memcmp(ampseq, seq, seqlen))
	      )
	    {
	      /* update info */
	      bp->swarmid = seed->swarmid;
	      bp->generation = seed->generation + 1;

	      /* lock mutex before adding this amplicon to swarm */
	      pthread_mutex_lock(&mutex_varmatch);

	      /* add to swarm */
	      current_swarm_tail->swarm_next = bp;
	      current_swarm_tail = bp;
	      
	      /* unlock mutex after adding this amplicon to swarm */
	      pthread_mutex_unlock(&mutex_varmatch);
	    }
	}
      j = (j + 1) % amphashtablesize;
    }
}

void generate_variants(unsigned long thread,
		       struct bucket_s * seed,
		       unsigned long start,
		       unsigned long len)
{
  /* 
     Generate all possible variants involving mutations from position start
     and extending len nucleotides. Insertions in front of those positions
     are included, but not those after. Positions are zero-based.
     The range may extend beyond the the length of the sequence indicating
     that inserts at the end of the sequence should be generated.

     The last thread will handle insertions at the end of the sequence,
     as well as identical sequences (no mutations).
  */

  unsigned char * varseq = ti[thread].varseq;
  
  unsigned char * seq = (unsigned char*) db_getsequence(seed->seqno);
  unsigned long seqlen = db_getsequencelen(seed->seqno);
  unsigned long end = MIN(seqlen,start+len);

  /* make an exact copy */
  memcpy(varseq, seq, seqlen);
  
#if 1
  /* identical non-variant */
  if (thread == threads -1)
    find_variant_matches(thread, varseq, seqlen, seed);
#endif

  /* substitutions */
  for(int i=start; i<end; i++)
    {
      for (int v=1; v<5; v++)
	if (v != seq[i])
	  {
	    varseq[i] = v;
	    find_variant_matches(thread, varseq, seqlen, seed);
	  }
      varseq[i] = seq[i];
    }

  /* deletions */
  memcpy(varseq, seq, start);
  if (start < seqlen-1)
    memcpy(varseq+start, seq+start+1, seqlen-start-1);
  for(int i=start; i<end; i++)
    {
      if ((i==0) || (seq[i] != seq[i-1]))
	find_variant_matches(thread, varseq, seqlen-1, seed);      
      varseq[i] = seq[i];
    }
  
  /* insertions */
  memcpy(varseq, seq, start);
  memcpy(varseq+start+1, seq+start, seqlen-start);
  for(int i=start; i<start+len; i++)
    {
      for(int v=1; v<5; v++)
	if((i==seqlen) || (v != seq[i]))
	  {
	    varseq[i] = v;
	    find_variant_matches(thread, varseq, seqlen+1, seed);
	  }
      if (i<seqlen)
	varseq[i] = seq[i];
    }
}

void * worker(void * vp)
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
	  generate_variants(t, tip->seed, tip->mut_start, tip->mut_length);
	  tip->work = 0;
	  pthread_cond_signal(&tip->workcond);
	}
    }

  pthread_mutex_unlock(&tip->workmutex);
  return 0;
}

void process_seed(struct bucket_s * seed)
{
  unsigned long seqlen = db_getsequencelen(seed->seqno);

  unsigned long thr = threads;
  if (thr > seqlen + 1)
    thr = seqlen+1;

  /* prepare work for the threads */
  unsigned long start = 0;
  for(unsigned long t=0; t<thr; t++)
    {
      struct thread_info_s * tip = ti + t;
      unsigned long length = (seqlen - start + thr - t) / (thr - t);
      tip->seed = seed;
      tip->mut_start = start;
      tip->mut_length = length;
      start += length;
      
      pthread_mutex_lock(&tip->workmutex);
      tip->work = 1;
      pthread_cond_signal(&tip->workcond);
      pthread_mutex_unlock(&tip->workmutex);
    }

  /* wait for theads to finish their work */
  for(int t=0; t<thr; t++)
    {
      struct thread_info_s * tip = ti + t;
      pthread_mutex_lock(&tip->workmutex);
      while (tip->work > 0)
	pthread_cond_wait(&tip->workcond, &tip->workmutex);
      pthread_mutex_unlock(&tip->workmutex);
    }
}

void threads_init()
{
  pthread_mutex_init(&mutex_varmatch, NULL);
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
	
  /* allocate memory for thread info, incl the variant sequences */
  unsigned long longestamplicon = db_getlongestsequence();
  ti = (struct thread_info_s *) xmalloc(threads * sizeof(struct thread_info_s));
  
  /* init and create worker threads */
  for(int t=0; t<threads; t++)
    {
      struct thread_info_s * tip = ti + t;
      tip->varseq = (unsigned char*) xmalloc(longestamplicon+1);
      tip->work = 0;
      pthread_mutex_init(&tip->workmutex, NULL);
      pthread_cond_init(&tip->workcond, NULL);
      if (pthread_create(&tip->pthread, &attr, worker, (void*)(long)t))
	fatal("Cannot create thread");
    }
}

void threads_done()
{
  /* finish and clean up worker threads */
  for(int t=0; t<threads; t++)
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
      free(tip->varseq);
    }

  free(ti);

  pthread_attr_destroy(&attr);
  pthread_mutex_destroy(&mutex_varmatch);
}

void update_stats(struct bucket_s * bp)
{
  /* update swarm stats */
  swarmsize++;
  if (bp->generation > swarm_maxgen)
    swarm_maxgen = bp->generation;
  unsigned long abundance = db_getabundance(bp->seqno);
  abundance_sum += abundance;
  if (abundance == 1)
    singletons++;
}

void algo_d1_run()
{
  unsigned long longestamplicon = db_getlongestsequence();
  unsigned long amplicons = db_getsequencecount();

  threads_init();

  /* compute hash for all amplicons and store them in hash table */
  amphashtablesize = HASHFILLFACTOR * amplicons;
  amphashtable = (struct bucket_s *) xmalloc (amphashtablesize * sizeof(struct bucket_s));
  for(unsigned long i=0; i<amphashtablesize; i++)
    {
      amphashtable[i].seqno = -1;
      amphashtable[i].swarmid = 0;
    }
  swarms_head = 0;
  swarms_tail = 0;
  all_head = 0;
  all_tail = 0;
  for(int i=0; i<amplicons; i++)
    {
      unsigned long seqlen = db_getsequencelen(i);
      unsigned char * seq = (unsigned char *) db_getsequence(i);
      unsigned long hash = HASH(seq, seqlen);
      unsigned long j = hash % amphashtablesize;
      
      /* find the first empty bucket */
      while (amphashtable[j].seqno >= 0)
	j = (j + 1) % amphashtablesize;

      struct bucket_s * bp = amphashtable + j;
      bp->seqno = i;
      bp->hash = hash;
      bp->swarmid = 0;
      bp->all_next = 0;
      bp->swarms_next = 0;
      bp->swarm_next = 0;

      if (i==0)
	all_head = bp;
      else
	all_tail->all_next = bp;
      all_tail = bp;
    }

  unsigned char * dir = 0;
  unsigned long * hearray = 0;

  if (uclustfile)
    {
      dir = (unsigned char *) xmalloc(longestamplicon*longestamplicon);
      hearray = (unsigned long *) xmalloc(2 * longestamplicon * sizeof(unsigned long));
    }
  
  /* for each non-swarmed amplicon look for subseeds ... */
  unsigned long swarmid = 0;
  for(seed = all_head; seed; seed = seed->all_next)
    {
      if (seed->swarmid == 0)
	{
	  /* start a new swarm with a new initial seed */
	  swarmid++;
	  seed->swarmid = swarmid;
	  seed->generation = 0;
	  seed->swarm_next = 0;
	  seed->swarms_next = 0;

	  /* link up this initial seed in the list of swarms */
	  if (swarmid == 1)
	    swarms_head = seed;
	  else
	    swarms_tail->swarms_next = seed;
	  swarms_tail = seed;
	  current_swarm_tail = seed;

	  /* initialize swarm stats */
	  swarmsize = 0;
	  swarm_maxgen = 0;
	  abundance_sum = 0;
	  singletons = 0;

	  update_stats(seed);
	  
	  /* find the first generation matches */
	  process_seed(seed);

	  /* find later generation matches */
	  struct bucket_s * subseed = seed->swarm_next;
	  while(subseed)
	    {
	      process_seed(subseed);
	      subseed = subseed->swarm_next;
	    }

	  /* update statistics */
	  for (struct bucket_s * bp = seed->swarm_next; 
	       bp;
	       bp = bp->swarm_next)
	    update_stats(bp);

	  /* update overall statistics */
	  if (swarmsize > largest)
 	    largest = swarmsize;
	  if (swarm_maxgen > maxgen)
	    maxgen = swarm_maxgen;

	  /* output statistics to file */
	  if (statsfile)
	    {
	      fprintf(statsfile, "%lu\t%lu\t", swarmsize, abundance_sum);
	      fprint_id_noabundance(statsfile, seed->seqno);
	      fprintf(statsfile, "\t%lu\t%lu\t%lu\t%lu\n", 
		      db_getabundance(seed->seqno),
		      singletons, swarm_maxgen, swarm_maxgen);
	    }

	  /* output results for one swarm in native format */
	  if (!mothur)
	    {
	      for (struct bucket_s * bp = seed; bp; bp = bp->swarm_next)
		{
		  if (bp != seed)
		    fputc(SEPCHAR, outfile);
		  fprint_id(outfile, bp->seqno);
		}
	      fputc('\n', outfile);
	    }
      
	  /* output break_swarms info */
	  if (break_swarms)
	    for (struct bucket_s * bp = seed->swarm_next; 
		 bp;
		 bp = bp->swarm_next)
	      {
		fprintf(stderr, "@@\t");
		fprint_id_noabundance(stderr, seed->seqno);
		fprintf(stderr, "\t");
		fprint_id_noabundance(stderr, bp->seqno);
		fprintf(stderr, "\t%d\n", 1);
	      }
	  
	  /* output swarm in uclust format */
	  if (uclustfile)
	    {
	      fprintf(uclustfile, "C\t%u\t%lu\t*\t*\t*\t*\t*\t",
		      seed->swarmid-1, swarmsize);
	      fprint_id(uclustfile, seed->seqno);
	      fprintf(uclustfile, "\t*\n");
          
	      fprintf(uclustfile, "S\t%u\t%lu\t*\t*\t*\t*\t*\t",
		      seed->swarmid-1, db_getsequencelen(seed->seqno));
	      fprint_id(uclustfile, seed->seqno);
	      fprintf(uclustfile, "\t*\n");

	      for (struct bucket_s * bp = seed->swarm_next; bp; bp = bp->swarm_next)
		{
		  char * dseq = db_getsequence(bp->seqno);
		  char * dend = dseq + db_getsequencelen(bp->seqno);
		  char * qseq = db_getsequence(seed->seqno);
		  char * qend = qseq + db_getsequencelen(seed->seqno);

		  unsigned long nwscore = 0;
		  unsigned long nwdiff = 0;
		  char * nwalignment = NULL;
		  unsigned long nwalignmentlength = 0;

		  nw(dseq, dend, qseq, qend,
		     score_matrix_63, gapopen, gapextend,
		     & nwscore, & nwdiff, & nwalignmentlength, & nwalignment,
		     dir, hearray, 0, 0);
              
		  double percentid = 100.0 * (nwalignmentlength - nwdiff) /
		    nwalignmentlength;
              
		  fprintf(uclustfile,
			  "H\t%u\t%lu\t%.1f\t+\t0\t0\t%s\t",
			  seed->swarmid-1,
			  db_getsequencelen(bp->seqno),
			  percentid, 
			  nwdiff > 0 ? nwalignment : "=");
              
		  fprint_id(uclustfile, bp->seqno);
		  fprintf(uclustfile, "\t");
		  fprint_id(uclustfile, seed->seqno);
		  fprintf(uclustfile, "\n");
		  
		  if (nwalignment)
		    free(nwalignment);
		}
	    }

	}
    }

  unsigned long swarmcount = swarmid;

  /* dump swarms in mothur format */
  /* cannot do it earlier because we need to know the number of swarms */
  if (mothur)
    {
      fprintf(outfile, "swarm_%ld\t%lu", resolution, swarmcount);
      for (struct bucket_s * seed = all_head; seed; seed = seed->swarms_next)
	for (struct bucket_s * bp = seed; bp; bp = bp->swarm_next)
	  {
	    if (bp == seed)
	      fputc('\t', outfile);
	    else
	      fputc(',', outfile);
	    fprint_id(outfile, bp->seqno);
	  }
      fputc('\n', outfile);
    }
  
  fprintf(stderr, "\n");
  fprintf(stderr, "Number of swarms:  %lu\n", swarmid);
  fprintf(stderr, "Largest swarm:     %lu\n", largest);
  fprintf(stderr, "Max generations:   %lu\n", maxgen);

  threads_done();

  free(amphashtable);

  if (uclustfile)
    {
      free(dir);
      free(hearray);
    }
}
