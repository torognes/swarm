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
  look for all possible variants of a sequence with only one difference
*/

#include "swarm.h"

//#define DEBUG

#define SEPCHAR ' '
#define HASH hash_djb2
//#define HASH hash_djb2a
//#define HASH hash_fnv_1a_32
//#define HASH hash_fnv_1a_64

#define HASHFILLFACTOR 2

struct bucket_s
{
  long seqno;
  unsigned long hash;
  unsigned long swarmid;
  unsigned long generation;
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
static unsigned char * varseq = 0;

/* overall statistics */
static unsigned long maxgen = 0;
static unsigned long largest = 0;

/* per swarm statistics */
static unsigned long singletons = 0;
static unsigned long abundance_sum = 0;
static unsigned long swarmsize = 0;
static unsigned long swarm_maxgen = 0;

#if 0
static unsigned long collisions = 0;
static unsigned long accesses = 0;
#endif

int amp_compare(const void * a, const void * b)
{
  struct bucket_s * x = * (struct bucket_s **) a;
  struct bucket_s * y = * (struct bucket_s **) b;

  if (x->seqno < y->seqno)
    return -1;
  else if (x->seqno > y->seqno)
    return +1;
  else
    return 0;
}


void dumpseq(unsigned long seqno)
{
  unsigned long seqlen = db_getsequencelen(seqno);
  char * seq = db_getsequence(seqno);
  for(int j=0; j<seqlen; j++)
    putchar(sym_nt[(unsigned int)seq[j]]);
  printf("\n");
}

void find_variant_matches(unsigned char * seq,
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
#if 0
	  accesses++;
#endif
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
	      
	      /* add to swarm */
	      current_swarm_tail->swarm_next = bp;
	      bp->swarm_next = 0;
	      current_swarm_tail = bp;
	      
	      /* update swarm stats */
	      swarmsize++;
	      if (bp->generation > swarm_maxgen)
		swarm_maxgen = bp->generation;
	      unsigned long abundance = db_getabundance(bp->seqno);
	      abundance_sum += abundance;
	      if (abundance == 1)
		singletons++;
	      
	      /* output break_swarms info */
	      if (break_swarms)
		{
		  fprintf(stderr, "@@\t");
		  fprint_id_noabundance(stderr, seed->seqno);
		  fprintf(stderr, "\t");
		  fprint_id_noabundance(stderr, bp->seqno);
		  fprintf(stderr, "\t1\n");
		}
	    }
	}
#if 0
      else
	collisions++;
#endif
      j = (j + 1) % amphashtablesize;
    }
}

/* 
   process_seed:
   - generate sequences with 0 or 1 difference from the given sequence
   with each variant:
   - find non-swarmed matching amplicons
   then for all hits of variants:
   - optionally, sort hits by decreasing abundance (ampliconindex)
   - add them to swarm and update amplicons
   - output swarm_breaker info
*/

void process_seed(struct bucket_s * seed)
{
  unsigned char * seq = (unsigned char*) db_getsequence(seed->seqno);
  unsigned long seqlen = db_getsequencelen(seed->seqno);
  memcpy(varseq, seq, seqlen);
  
#if 1
  /* identical */
  find_variant_matches(varseq, seqlen, seed);
#endif

  /* substitutions */
  for(int i=0; i<seqlen; i++)
    {
      for (int v=1; v<5; v++)
	if (v != seq[i])
	  {
	    varseq[i] = v;
	    find_variant_matches(varseq, seqlen, seed);
	  }
      varseq[i] = seq[i];
    }

  /* deletions */
  for(int i=0; i<seqlen; i++)
    {
      if(i>0)
	{
	  if (varseq[i] != seq[i-1])
	    varseq[i] = seq[i-1];
	  else
	    continue;
	}
      find_variant_matches(varseq+1, seqlen-1, seed);
    }

  /* insertions */
  varseq[seqlen] = seq[seqlen-1];
  for(int i=0; i<=seqlen; i++)
    {
      for(int v=1; v<5; v++)
	{
	  if((i==seqlen) || (v != seq[i]))
	    {
	      varseq[i] = v;
	      find_variant_matches(varseq, seqlen+1, seed);
	    }
	}
      if (i<seqlen)
	varseq[i] = seq[i];
    }
}


void algo_d1_run()
{
  unsigned long longestamplicon = db_getlongestsequence();
  unsigned long amplicons = db_getsequencecount();
 
  varseq = (unsigned char*) xmalloc(longestamplicon+1);

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
	  swarmsize = 1;
	  swarm_maxgen = 0;
	  unsigned long seedabundance = db_getabundance(seed->seqno);
	  abundance_sum = seedabundance;
	  if (seedabundance == 1)
	    singletons = 1;
	  else
	    singletons = 0;
	  
	  /* find the first generation matches */
	  process_seed(seed);

	  /* find later generation matches */
	  struct bucket_s * subseed = seed->swarm_next;
	  while(subseed)
	    {
	      process_seed(subseed);
	      subseed = subseed->swarm_next;
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
      
	  /* update overall statistics */
	  if (swarmsize > largest)
 	    largest = swarmsize;
	  if (swarm_maxgen > maxgen)
	    maxgen = swarm_maxgen;

	  /* output swarm in uclust format */
	  if (uclustfile)
	    {
	      fprintf(uclustfile, "C\t%lu\t%lu\t*\t*\t*\t*\t*\t",
		      seed->swarmid-1, swarmsize);
	      fprint_id(uclustfile, seed->seqno);
	      fprintf(uclustfile, "\t*\n");
          
	      fprintf(uclustfile, "S\t%lu\t%lu\t*\t*\t*\t*\t*\t",
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
			  "H\t%lu\t%lu\t%.1f\t+\t0\t0\t%s\t",
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

	  if (statsfile)
	    {
	      fprintf(statsfile, "%lu\t%lu\t", swarmsize, abundance_sum);
	      fprint_id_noabundance(statsfile, seed->seqno);
	      fprintf(statsfile, "\t%lu\t%lu\t%lu\t%lu\n", 
		      seedabundance, singletons, swarm_maxgen, swarm_maxgen);
	    }
	}
    }

  unsigned long swarmcount = swarmid;

  /* dump swarms in mothur format */
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

#if 0
  fprintf(stderr, "\nHash collisions:   %lu\n", collisions);
  fprintf(stderr, "Hash accesses:     %lu\n", accesses);
#endif

  free(varseq);
  free(amphashtable);

  if (uclustfile)
    {
      free(dir);
      free(hearray);
    }
}
