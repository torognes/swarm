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

#define HASH hash_cityhash64

//#define REVCOMP

struct bucket
{
  uint64_t hash;
  unsigned int seqno_first;
  unsigned int seqno_last;
  uint64_t mass;
  unsigned int size;
  unsigned int singletons;
};

int derep_compare(const void * a, const void * b);

int derep_compare(const void * a, const void * b)
{
  const struct bucket * x = static_cast<const struct bucket *>(a);
  const struct bucket * y = static_cast<const struct bucket *>(b);

  /* highest abundance first, otherwise keep order */

  if (x->mass < y->mass)
    return +1;
  else if (x->mass > y->mass)
    return -1;
  else
    {
      if (x->seqno_first < y->seqno_first)
        return -1;
      else if (x->seqno_first > y->seqno_first)
        return +1;
      else
        return 0;
    }
}

#ifdef REVCOMP
char map_complement[5] = { 0, 4, 3, 2, 1 };

void reverse_complement(char * rc, char * seq, int64_t len)
{
  /* Write the reverse complementary sequence to rc.
     The memory for rc must be long enough for the rc of the sequence
     (identical to the length of seq + 1). */

  for(int64_t i=0; i<len; i++)
    rc[i] = map_complement[(int)(1 + nt_extract(seq, len-1-i))];
  rc[len] = 0;
}
#endif

void dereplicate()
{
  /* adjust size of hash table for 2/3 fill rate */
  uint64_t dbsequencecount = db_getsequencecount();
  uint64_t hashtablesize = 1;
  while (100 * dbsequencecount > 70 * hashtablesize)
    hashtablesize <<= 1;
  uint64_t derep_hash_mask = hashtablesize - 1;

  struct bucket * hashtable =
    static_cast<struct bucket *>(xmalloc(sizeof(bucket) * hashtablesize));

  memset(hashtable, 0, sizeof(bucket) * hashtablesize);

  uint64_t swarmcount = 0;
  uint64_t maxmass = 0;
  unsigned int maxsize = 0;

  /* alloc and init table of links to other sequences in cluster */
  unsigned int * nextseqtab = static_cast<unsigned int *>
    (xmalloc(sizeof(unsigned int) * dbsequencecount));
  memset(nextseqtab, 0, sizeof(unsigned int) * dbsequencecount);

#ifdef REVCOMP
  /* allocate memory for reverse complementary sequence */
  char * rc_seq = (char*) xmalloc(db_getlongestsequence() + 1);
#endif

  progress_init("Dereplicating:    ", dbsequencecount);

  for(unsigned int i=0; i<dbsequencecount; i++)
    {
      unsigned int seqlen = db_getsequencelen(i);
      char * seq = db_getsequence(i);

      /*
        Find free bucket or bucket for identical sequence.
        Make sure sequences are exactly identical
        in case of any hash collision.
        With 64-bit hashes, there is about 50% chance of a
        collision when the number of sequences is about 5e9.
      */

      uint64_t hash = HASH(reinterpret_cast<unsigned char *>(seq),
                           nt_bytelength(seqlen));
      uint64_t j = hash & derep_hash_mask;
      struct bucket * bp = hashtable + j;

      while ((bp->mass) &&
             ((bp->hash != hash) ||
              (seqlen != db_getsequencelen(bp->seqno_first)) ||
              (memcmp(seq,
                      db_getsequence(bp->seqno_first),
                      nt_bytelength(seqlen)))))
        {
          bp++;
          j++;
          if (bp >= hashtable + hashtablesize)
            {
              bp = hashtable;
              j = 0;
            }
        }

#ifdef REVCOMP
      if (! bp->mass)
        {
          /* no match on plus strand */
          /* check minus strand as well */

          reverse_complement(rc_seq, seq, seqlen);
          uint64_t rc_hash = HASH((unsigned char*)rc_seq, nt_bytelength(seqlen));
          struct bucket * rc_bp = hashtable + rc_hash % hashtablesize;
          uint64_t k = rc_hash & derep_hash_mask;

          while ((rc_bp->mass) &&
                 ((rc_bp->hash != rc_hash) ||
                  (seqlen != db_getsequencelen(rc_bp->seqno_first)) ||
                  (memcmp(rc_seq,
                          db_getsequence(rc_bp->seqno_first),
                          nt_bytelength(seqlen)))))
            {
              rc_bp++;
              k++;
              if (rc_bp >= hashtable + hashtablesize)
                {
                  rc_bp = hashtable;
                  k++;
                }
            }

          if (rc_bp->mass)
            {
              bp = rc_bp;
              j = k;
            }
        }
#endif

      uint64_t ab = db_getabundance(i);

      if (bp->mass)
        {
          /* at least one identical sequence already */
          nextseqtab[bp->seqno_last] = i;
        }
      else
        {
          /* no identical sequences yet, start a new cluster */
          swarmcount++;
          bp->hash = hash;
          bp->seqno_first = i;
          bp->size = 0;
          bp->singletons = 0;
        }

      bp->size++;
      bp->seqno_last = i;
      bp->mass += ab;

      if (ab == 1)
        bp->singletons++;

      if (bp->mass > maxmass)
        maxmass = bp->mass;

      if (bp->size > maxsize)
        maxsize = bp->size;

      progress_update(i);
    }
  progress_done();

#ifdef REVCOMP
  xfree(rc_seq);
#endif

  progress_init("Sorting:          ", 1);
  qsort(hashtable, hashtablesize, sizeof(bucket), derep_compare);
  progress_done();


  /* dump swarms */

  progress_init("Writing swarms:   ", swarmcount);

  if (opt_mothur)
    fprintf(outfile, "swarm_%" PRId64 "\t%" PRIu64, opt_differences, swarmcount);

  for(unsigned int i = 0; i < swarmcount; i++)
    {
      unsigned int seed = hashtable[i].seqno_first;
      if (opt_mothur)
        fputc('\t', outfile);
      fprint_id(outfile, seed);
      unsigned int a = nextseqtab[seed];

      while (a)
        {
          if (opt_mothur)
            fputc(',', outfile);
          else
            fputc(SEPCHAR, outfile);
          fprint_id(outfile, a);
          a = nextseqtab[a];
        }

      if (!opt_mothur)
        fputc('\n', outfile);

      progress_update(i+1);
    }

  if (opt_mothur)
    fputc('\n', outfile);

  progress_done();


  /* dump seeds in fasta format with sum of abundances */

  if (opt_seeds)
    {
      progress_init("Writing seeds:    ", swarmcount);
      for(unsigned int i=0; i < swarmcount; i++)
        {
          unsigned int seed = hashtable[i].seqno_first;
          fprintf(fp_seeds, ">");
          fprint_id_with_new_abundance(fp_seeds, seed, hashtable[i].mass);
          fprintf(fp_seeds, "\n");
          db_fprintseq(fp_seeds, seed, 0);
          progress_update(i+1);
        }
      progress_done();
    }

  /* output swarm in uclust format */

  if (uclustfile)
    {
      progress_init("Writing UCLUST:   ", swarmcount);

      for(unsigned int swarmid = 0; swarmid < swarmcount ; swarmid++)
        {
          struct bucket * bp = hashtable + swarmid;

          unsigned int seed = bp->seqno_first;

          fprintf(uclustfile, "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                  swarmid,
                  bp->size);
          fprint_id(uclustfile, seed);
          fprintf(uclustfile, "\t*\n");

          fprintf(uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                  swarmid,
                  db_getsequencelen(seed));
          fprint_id(uclustfile, seed);
          fprintf(uclustfile, "\t*\n");

          unsigned int a = nextseqtab[seed];

          while (a)
            {
              fprintf(uclustfile,
                      "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                      swarmid,
                      db_getsequencelen(a),
                      100.0,
                      "=");
              fprint_id(uclustfile, a);
              fprintf(uclustfile, "\t");
              fprint_id(uclustfile, seed);
              fprintf(uclustfile, "\n");
              a = nextseqtab[a];
            }

          progress_update(swarmid+1);
        }
      progress_done();
    }

  /* output internal structure to file */

  if (opt_internal_structure)
    {
      progress_init("Writing structure:", swarmcount);

      for(uint64_t i = 0; i < swarmcount; i++)
        {
          struct bucket * sp = hashtable + i;
          uint64_t seed = sp->seqno_first;
          unsigned int a = nextseqtab[seed];
          while (a)
            {
              fprint_id_noabundance(internal_structure_file, seed);
              fprintf(internal_structure_file, "\t");
              fprint_id_noabundance(internal_structure_file, a);
              fprintf(internal_structure_file, "\t%d\t%" PRIu64 "\t%d\n", 0, i+1, 0);
              a = nextseqtab[a];
            }
          progress_update(i);
        }
      progress_done();
    }

  /* output statistics to file */

  if (statsfile)
    {
      progress_init("Writing stats:    ", swarmcount);
      for(uint64_t i = 0; i < swarmcount; i++)
        {
          struct bucket * sp = hashtable + i;
          fprintf(statsfile, "%u\t%" PRIu64 "\t", sp->size, sp->mass);
          fprint_id_noabundance(statsfile, sp->seqno_first);
          fprintf(statsfile, "\t%" PRIu64 "\t%u\t%u\t%u\n",
                  db_getabundance(sp->seqno_first),
                  sp->singletons, 0U, 0U);
          progress_update(i);
        }
      progress_done();
    }


  fprintf(logfile, "\n");
  fprintf(logfile, "Number of swarms:  %" PRIu64 "\n", swarmcount);
  fprintf(logfile, "Largest swarm:     %u\n", maxsize);
  fprintf(logfile, "Heaviest swarm:    %" PRIu64 "\n", maxmass);

  xfree(nextseqtab);
  xfree(hashtable);
}
