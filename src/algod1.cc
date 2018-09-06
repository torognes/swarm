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

/*
  This version of the swarm algorithm uses Frederic's idea for d=1 to
  enumerate all of the maximum 7L+4 possible variants of a sequence
  with only one difference, where L is the length of the sequence.
*/

#include "swarm.h"

#define HASHFILLFACTOR 0.5
//#define HASHSTATS

/* Information about each amplicon */

static struct ampinfo_s
{
  int swarmid;
  int parent;
  int generation;
  int next;       /* amp id of next amplicon in swarm */
  int graft_cand; /* amp id of potential grafting parent (fastidious) */
} * ampinfo = 0;

/* Information about each swarm (OTU) */

static struct swarminfo_s
{
  int seed; /* amplicon id of the initial seed of this swarm */
  int last; /* amplicon id of the last seed in this swarm */
  int size; /* total number of amplicons in this swarm */
  int singletons; /* number of amplicons with abundance 1 */
  int maxgen; /* the generation of the amplicon farthest from the seed */
  long mass; /* the sum of abundances of amplicons in this swarm */
  long sumlen; /* sum of length of amplicons in swarm */
  bool attached; /* this is a small swarm attached to a large (fastidious) */
} * swarminfo = 0;

static long swarminfo_alloc = 0;

/* Information about potential grafts */
static long graft_candidates = 0;
static pthread_mutex_t graft_mutex;

#define NO_SWARM (-1)

static int current_swarm_tail;

static unsigned long hash_tablesize = 0;

/* overall statistics */
static int maxgen = 0;
static int largest = 0;
static long swarmcount_adjusted = 0;

/* per swarm statistics */
static unsigned long singletons = 0;
static unsigned long abundance_sum = 0; /* = mass */
static int swarmsize = 0;
static int swarm_maxgen = 0;
static unsigned long swarm_sumlen = 0;

static struct thread_info_s
{
  pthread_t pthread;
  pthread_mutex_t workmutex;
  pthread_cond_t workcond;
  int work;
  unsigned char * varseq;
  int seed;
  unsigned long mut_start;
  unsigned long mut_length;
  int * hits_data;
  int hits_alloc;
  int hits_count;
} * ti;

static pthread_attr_t attr;

#ifdef HASHSTATS
unsigned long hits = 0;
unsigned long bloom_matches = 0;
unsigned long success = 0;
unsigned long tries  = 0;
unsigned long bingo = 0;
unsigned long collisions = 0;
#endif

static int hash_shift;
static unsigned long hash_mask;
static unsigned char * hash_occupied = 0;
static unsigned long * hash_values = 0;
static int * hash_data = 0;

#define BLOOM_PATTERN_BITS 10
#define BLOOM_PATTERN_COUNT (1 << BLOOM_PATTERN_BITS)

static unsigned long * hash_bloom = 0;
static unsigned long bit_patterns[BLOOM_PATTERN_COUNT];
static unsigned long bloom_size_bytes = 0;
static unsigned long bloom_mask = 0;
static unsigned long bloom_pattern_mask = 0;

static int * global_hits_data = 0;
static int global_hits_alloc = 0;
static int global_hits_count = 0;

static unsigned long threads_used = 0;

void generate_bit_patterns()
{
  for (unsigned int i = 0; i < BLOOM_PATTERN_COUNT; i++)
    {
      unsigned long pattern = 0;
      for (unsigned int j = 0; j < 8; j++)
        {
          unsigned long onebit = 1ULL << (random() & 63);
          while (pattern & onebit)
            onebit = 1ULL << (random() & 63);
          pattern |= onebit;
        }
      bit_patterns[i] = pattern;
    }
}

void bloom_init()
{
  bloom_size_bytes = hash_tablesize;
  bloom_mask = (bloom_size_bytes >> 3) - 1;
  hash_bloom = (unsigned long *) xmalloc(bloom_size_bytes);
  memset(hash_bloom, 0xff, bloom_size_bytes);

  bloom_pattern_mask = BLOOM_PATTERN_COUNT - 1;
  generate_bit_patterns();
}

void bloom_exit()
{
  free(hash_bloom);
}

inline void bloom_set(unsigned long hash)
{
  unsigned long bloom_pattern = bit_patterns[hash & bloom_pattern_mask];
  hash_bloom[(hash >> BLOOM_PATTERN_BITS) & bloom_mask] &= ~bloom_pattern;
}

inline bool bloom_get(unsigned long hash)
{
  unsigned long bloom_pattern = bit_patterns[hash & bloom_pattern_mask];
  unsigned long bloom_bits
    = hash_bloom[(hash >> BLOOM_PATTERN_BITS) & bloom_mask];
  return (bloom_bits & bloom_pattern) == 0;
}

inline unsigned int hash_getindex(unsigned long hash)
{
  // Shift bits right to get independence from the simple Bloom filter hash
  hash = hash >> 38;
  return hash & hash_mask;
}

inline unsigned int hash_getnextindex(unsigned int j)
{
  return (j+1) & hash_mask;
}

void hash_alloc(unsigned long amplicons)
{
  hash_tablesize = 1;
  hash_shift = 0;
  while (amplicons > HASHFILLFACTOR * hash_tablesize)
    {
      hash_tablesize <<= 1;
      hash_shift++;
    }
  hash_mask = hash_tablesize - 1;

  hash_occupied =
    (unsigned char *) xmalloc((hash_tablesize + 63) / 8);
  memset(hash_occupied, 0, (hash_tablesize + 63) / 8);

  hash_values =
    (unsigned long *) xmalloc(hash_tablesize * sizeof(unsigned long));

  hash_data =
    (int *) xmalloc(hash_tablesize * sizeof(int));

  bloom_init();
}


void hash_free()
{
  bloom_exit();

  free(hash_occupied);
  free(hash_values);
  free(hash_data);
}

inline void hash_set_occupied(unsigned int j)
{
  hash_occupied[j >> 3] |= (1 << (j & 7));
}

inline int hash_is_occupied(unsigned int j)
{
  return hash_occupied[j >> 3] & (1 << (j & 7));
}

inline void hash_set_value(unsigned int j, unsigned long hash)
{
  hash_values[j] = hash;
}

inline int hash_compare_value(unsigned int j, unsigned long hash)
{
  return (hash_values[j] == hash);
}

inline void hash_insert(int amp)
{
  /* find the first empty bucket */
  unsigned long hash = db_gethash(amp);
  unsigned int j = hash_getindex(hash);
  while (hash_is_occupied(j))
    j = hash_getnextindex(j);

  hash_set_occupied(j);
  hash_set_value(j, hash);
  hash_data[j] = amp;

  bloom_set(hash);
}

enum variant_enum
  {
    original,
    substitution,
    deletion,
    insertion
  };

bool seq_identical(char * a,
                   int a_start,
                   char * b,
                   int b_start,
                   int length)
{
  /* compare parts of two compressed sequences a and b */
  /* return false if different, true if identical */

  for(int i = 0; i < length; i++)
    if (nt_extract(a, a_start + i) != nt_extract(b, b_start + i))
      return false;
  return true;
}

inline bool check_variant(int seed,
                          variant_enum vartype,
                          unsigned long pos,
                          unsigned char base,
                          int amp)
{
  /* check if amp is identical to seed with the given variant */

  char * seed_sequence = db_getsequence(seed);
  unsigned int seed_seqlen = db_getsequencelen(seed);

  char * amp_sequence = db_getsequence(amp);
  unsigned int amp_seqlen = db_getsequencelen(amp);

  switch (vartype)
    {

    case original:

      if (seed_seqlen != amp_seqlen)
        return false;

      if (! seq_identical(seed_sequence, 0,
                          amp_sequence, 0,
                          seed_seqlen))
        return false;

      break;

    case substitution:

      if (seed_seqlen != amp_seqlen)
        return false;

      if (nt_extract(amp_sequence, pos) != base)
        return false;

      if (! seq_identical(seed_sequence, 0,
                          amp_sequence, 0,
                          pos))
        return false;

      if (! seq_identical(seed_sequence, pos + 1,
                          amp_sequence,  pos + 1,
                          seed_seqlen - pos - 1))
        return false;

      break;

    case deletion:

      if ((seed_seqlen - 1) != amp_seqlen)
        return false;

      if (! seq_identical(seed_sequence, 0,
                          amp_sequence, 0,
                          pos))
        return false;

      if (! seq_identical(seed_sequence, pos + 1,
                          amp_sequence,  pos,
                          seed_seqlen - pos - 1))
        return false;
      break;

    case insertion:

      if ((seed_seqlen + 1) != amp_seqlen)
        return false;

      if (nt_extract(amp_sequence, pos) != base)
        return false;

      if (! seq_identical(seed_sequence, 0,
                          amp_sequence, 0,
                          pos))
        return false;

      if (! seq_identical(seed_sequence, pos,
                          amp_sequence,  pos + 1,
                          seed_seqlen - pos))
        return false;

      break;

    default:

      break;
    }

  return true;
}

void find_variant_matches2(unsigned long thread,
                           int seed,
                           unsigned long hash,
                           variant_enum vartype,
                           unsigned int pos,
                           unsigned int base,
                           unsigned long max_abundance)
{
#ifdef HASHSTATS
  bloom_matches++;
#endif

  /* compute hash and corresponding hash table index */

  unsigned int j = hash_getindex(hash);

  /* find matching buckets */

  struct thread_info_s * tip = ti + thread;

  while (hash_is_occupied(j))
    {
#ifdef HASHSTATS
      hits++;
#endif
      if (hash_compare_value(j, hash))
        {
#ifdef HASHSTATS
          success++;
#endif

          /* check if not already swarmed */
          int amp = hash_data[j];
          if ((ampinfo[amp].swarmid == NO_SWARM) &&
              (db_getabundance(amp) <= max_abundance))
            {
              if (check_variant(seed, vartype, pos, base, amp))
                {
#ifdef HASHSTATS
                  bingo++;
#endif
                  tip->hits_data[tip->hits_count++] = amp;
                  break;
                }
#ifdef HASHSTATS
              else
                {
                  collisions++;

#ifdef HASHDETAILS
                  fprintf(logfile, "Hash collision between ");
                  fprint_id_noabundance(logfile, seed);
                  fprintf(logfile, " and ");
                  fprint_id_noabundance(logfile, amp);
                  fprintf(logfile, ".\n");
#endif
                }
#endif
            }
        }
      j = hash_getnextindex(j);
    }
}

inline void find_variant_matches(unsigned long thread,
                                 int seed,
                                 unsigned long hash,
                                 variant_enum vartype,
                                 unsigned int pos,
                                 unsigned int base,
                                 unsigned long max_abundance)
{
#ifdef HASHSTATS
  tries++;
#endif

  if (bloom_get(hash))
    find_variant_matches2(thread,
                          seed,
                          hash,
                          vartype,
                          pos,
                          base,
                          max_abundance);
}

void generate_variants(unsigned long thread,
                       int seed,
                       unsigned long start,
                       unsigned long len)
{
  (void) start;
  (void) len;
  if (thread > 0)
    return;

  /*
     Generate all possible variants involving mutations from position start
     and extending len nucleotides. Insertions in front of those positions
     are included, but not those after. Positions are zero-based.
     The range may extend beyond the length of the sequence indicating
     that inserts at the end of the sequence should be generated.

     The last thread will handle insertions at the end of the sequence,
     as well as identical sequences (no mutations).
  */

  /* For now, thread 0 will do all work. Start and Len are ignored. */

  ti[thread].hits_count = 0;
  unsigned long m = opt_no_otu_breaking ? ULONG_MAX : db_getabundance(seed);
  char * sequence = db_getsequence(seed);
  unsigned int seqlen = db_getsequencelen(seed);
  unsigned long hash = db_gethash(seed);

  /* identical non-variant */

  find_variant_matches(thread, seed, hash, original, 0, 0, m);

  /* substitutions */

  for(unsigned int i = 0; i < seqlen; i++)
    {
      unsigned int base = nt_extract(sequence, i);
      unsigned long hash1 = hash ^ zobrist_value(i, base);
      for (unsigned int v = 0; v < 4; v ++)
        if (v != base)
          {
            unsigned long hash2 = hash1 ^ zobrist_value(i, v);
            find_variant_matches(thread, seed, hash2, substitution, i, v, m);
          }
    }

  /* deletions */

  hash = zobrist_hash_delete_first((unsigned char *) sequence, seqlen);
  find_variant_matches(thread, seed, hash, deletion, 0, 0, m);
  unsigned int base = nt_extract(sequence, 0);
  for(unsigned int i = 1; i < seqlen; i++)
    {
      unsigned int v = nt_extract(sequence, i);
      if (v != base)
        {
          hash ^= zobrist_value(i - 1, base) ^ zobrist_value(i - 1, v);
          find_variant_matches(thread, seed, hash, deletion, i, 0, m);
          base = v;
        }
    }

  /* insertions */

  hash = zobrist_hash_insert_first((unsigned char *) sequence, seqlen);
  for (unsigned int v = 0; v < 4; v++)
    {
      unsigned long hash1 = hash ^ zobrist_value(0, v);
      find_variant_matches(thread, seed, hash1, insertion, 0, v, m);
    }
  for (unsigned int i = 0; i < seqlen; i++)
    {
      unsigned int base = nt_extract(sequence, i);
      hash ^= zobrist_value(i, base) ^ zobrist_value(i+1, base);
      for (unsigned int v = 0; v < 4; v++)
        if (v != base)
          {
            unsigned long hash1 = hash ^ zobrist_value(i + 1, v);
            find_variant_matches(thread, seed, hash1, insertion, i + 1, v, m);
          }
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
      while (tip->work == 0)
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

void threads_init()
{
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* allocate memory for thread info, incl the variant sequences */
  unsigned long longestamplicon = db_getlongestsequence();
  ti = (struct thread_info_s *)
    xmalloc(opt_threads * sizeof(struct thread_info_s));

  /* init and create worker threads */
  for(long t=0; t<opt_threads; t++)
    {
      struct thread_info_s * tip = ti + t;
      tip->varseq = (unsigned char*) xmalloc(longestamplicon+1);
      tip->hits_alloc = 7 * longestamplicon + 4;
      tip->hits_data = (int*) xmalloc(tip->hits_alloc * sizeof(int));
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
      free(tip->varseq);
      free(tip->hits_data);
    }

  free(ti);

  pthread_attr_destroy(&attr);
}

void add_amp_to_swarm(int amp)
{
  /* add to swarm */
  ampinfo[current_swarm_tail].next = amp;
  current_swarm_tail = amp;
}

void process_seed(int subseed)
{
  unsigned long seqlen = db_getsequencelen(subseed);

  threads_used = opt_threads;
  if (threads_used > seqlen + 1)
    threads_used = seqlen+1;

  /* prepare work for the threads */
  unsigned long start = 0;
  for(unsigned long t=0; t<threads_used; t++)
    {
      struct thread_info_s * tip = ti + t;
      unsigned long length =
        (seqlen - start + threads_used - t) / (threads_used - t);
      tip->seed = subseed;
      tip->mut_start = start;
      tip->mut_length = length;
      start += length;

      pthread_mutex_lock(&tip->workmutex);
      tip->work = 1;
      pthread_cond_signal(&tip->workcond);
      pthread_mutex_unlock(&tip->workmutex);
    }

  /* wait for threads to finish their work */
  for(unsigned int t=0; t<threads_used; t++)
    {
      struct thread_info_s * tip = ti + t;
      pthread_mutex_lock(&tip->workmutex);
      while (tip->work > 0)
        pthread_cond_wait(&tip->workcond, &tip->workmutex);
      pthread_mutex_unlock(&tip->workmutex);
    }

  /* join hits from the threads */

  for(unsigned int t=0; t<threads_used; t++)
    {
      if (global_hits_count + ti[t].hits_count > global_hits_alloc)
        {
          while (global_hits_count + ti[t].hits_count > global_hits_alloc)
            global_hits_alloc <<= 1;
          global_hits_data = (int*)xrealloc(global_hits_data,
                                            global_hits_alloc * sizeof(int));
        }
      for(int i=0; i < ti[t].hits_count; i++)
        {
          long amp = ti[t].hits_data[i];

          /* add to list for this generation */
          global_hits_data[global_hits_count++] = amp;

          /* update info */
          ampinfo[amp].swarmid = ampinfo[subseed].swarmid;
          ampinfo[amp].generation = ampinfo[subseed].generation + 1;
          ampinfo[amp].parent = subseed;
        }
    }
}

void update_stats(int amp)
{
  /* update swarm stats */
  struct ampinfo_s * bp = ampinfo + amp;

  swarmsize++;
  if (bp->generation > swarm_maxgen)
    swarm_maxgen = bp->generation;
  unsigned long abundance = db_getabundance(amp);
  abundance_sum += abundance;
  if (abundance == 1)
    singletons++;
  swarm_sumlen += db_getsequencelen(amp);
}

void attach(int seed, int amp)
{
  /* graft light swarm (amp) on heavy swarm (seed) */

#if 0
  fprintf(logfile,
          "\nGrafting light swarm with amplicon %d on "
          "heavy swarm with amplicon %d (swarm ids: %d %d)\n",
          amp,
          seed,
          ampinfo[amp].swarmid,
          ampinfo[seed].swarmid);
#endif

  swarminfo_s * hp = swarminfo + ampinfo[seed].swarmid;
  swarminfo_s * lp = swarminfo + ampinfo[amp].swarmid;

  // attach the seed of the light swarm to the tail of the heavy swarm
  ampinfo[hp->last].next = lp->seed;
  hp->last = lp->last;

  // Update swarm info
  hp->size += lp->size;
  hp->singletons += lp->singletons;
  hp->mass += lp->mass;
  hp->sumlen += lp->sumlen;
  /* maxgen is untouched */

  /* flag attachment to avoid doing it again */
  lp->attached = true;

  // Update overall stats
  if (hp->size > largest)
    largest = hp->size;

  swarmcount_adjusted--;
}

void add_graft_candidate(int seed, int amp)
{
  pthread_mutex_lock(&graft_mutex);
  graft_candidates++;
  if ((ampinfo[amp].graft_cand == NO_SWARM)||(ampinfo[amp].graft_cand > seed))
    ampinfo[amp].graft_cand = seed;
  pthread_mutex_unlock(&graft_mutex);
}

struct graft_cand
{
  int parent;
  int child;
} * graft_array;

int compare_grafts(const void * a, const void * b)
{
  struct graft_cand * x = (struct graft_cand *) a;
  struct graft_cand * y = (struct graft_cand *) b;
  if (x->parent < y->parent)
    return -1;
  else if (x->parent > y->parent)
    return +1;
  else
    if (x->child < y->child)
      return -1;
    else if (x->child > y->child)
      return +1;
    else
      return 0;
}

int attach_candidates(int amplicons)
{
  /* count pairs */
  int pair_count = 0;
  for(int i=0; i < amplicons; i++)
    if (ampinfo[i].graft_cand != NO_SWARM)
      pair_count++;

  int grafts = 0;
  progress_init("Grafting light swarms on heavy swarms", pair_count);

  /* allocate memory */
  graft_array = (struct graft_cand *)
    xmalloc(pair_count * sizeof(struct graft_cand));

  /* fill in */
  int j = 0;
  for(int i=0; i < amplicons; i++)
    if (ampinfo[i].graft_cand != NO_SWARM)
      {
        graft_array[j].parent = ampinfo[i].graft_cand;
        graft_array[j].child = i;
        j++;
      }

  /* sort */
  qsort(graft_array, pair_count, sizeof(struct graft_cand), compare_grafts);

  /* attach in order */
  for(int i=0; i < pair_count; i++)
    {
      int parent = graft_array[i].parent;
      int child  = graft_array[i].child;

      if (swarminfo[ampinfo[child].swarmid].attached)
        {
          /* this light swarm is already attached */
          ampinfo[child].graft_cand = NO_SWARM;
        }
      else
        {
          /* attach child to parent */
          attach(parent, child);
          grafts++;
        }
      progress_update(i+1);
    }
  progress_done();
  free(graft_array);
  return grafts;
}

bool hash_check_attach(char * seq,
                       unsigned long seqlen,
                       int seed)
{
  /* compute hash and corresponding hash table index */

  unsigned long hash = zobrist_hash((unsigned char*)seq, seqlen);
  unsigned int j = hash_getindex(hash);

  /* find matching buckets */

  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, hash))
        {
          /* check that mass is below threshold */
          int amp = hash_data[j];

          struct swarminfo_s * smallp = swarminfo + ampinfo[amp].swarmid;

          if (smallp->mass < opt_boundary)
            {
              unsigned long ampseqlen = db_getsequencelen(amp);
              unsigned char * ampseq = (unsigned char *) db_getsequence(amp);

              /* make absolutely sure sequences are identical */
              if ((ampseqlen == seqlen) && (!memcmp(ampseq, seq, seqlen)))
                {
                  add_graft_candidate(seed, amp);
                  return 1;
                }
            }
        }
      j = hash_getnextindex(j);
    }
  return 0;
}

long expected_variant_count(char * seq, int len)
{
  int c = 0;
  for(int i=1; i<len; i++)
    if (seq[i] != seq[i-1])
      c++;
  return 6*len+5+c;
}


long fastidious_mark_small_var(BloomFilter * bloom,
                               char * varseq,
                               int seed)
{
  /*
    add all microvariants of seed to Bloom filter

    bloom is a BloomFilter in which to enter the variants
    buffer is a buffer large enough to hold all sequences + 1 insertion
    seed is the original seed
  */

  long variants = 0;

  unsigned char * seq = (unsigned char*) db_getsequence(seed);
  unsigned long seqlen = db_getsequencelen(seed);

  /* make an exact copy */
  memcpy(varseq, seq, seqlen);

  /* substitutions */
  for(unsigned int i=0; i<seqlen; i++)
    {
      for (int v=1; v<5; v++)
        if (v != seq[i])
          {
            varseq[i] = v;
            bloom->set(varseq, seqlen);
            variants++;
          }
      varseq[i] = seq[i];
    }

  /* deletions */
  if (seqlen > 1)
    memcpy(varseq, seq+1, seqlen-1);
  for(unsigned int i=0; i<seqlen; i++)
    {
      if ((i==0) || (seq[i] != seq[i-1]))
        {
          bloom->set(varseq, seqlen-1);
          variants++;
        }
      varseq[i] = seq[i];
    }

  /* insertions */
  memcpy(varseq+1, seq, seqlen);
  for(unsigned int i=0; i<seqlen+1; i++)
    {
      for(int v=1; v<5; v++)
        {
          if((i==seqlen) || (v != seq[i]))
            {
              varseq[i] = v;
              bloom->set(varseq, seqlen+1);
              variants++;
            }
        }
      if (i<seqlen)
        varseq[i] = seq[i];
    }
#if 0
  long e = expected_variant_count((char*)seq, seqlen);
  if (variants != e)
    fprintf(logfile, "Incorrect number of variants: %ld Expected: %ld\n", variants, e);
#endif
  return variants;
}

long fastidious_check_large_var_2(char * seq,
                                  size_t seqlen,
                                  char * varseq,
                                  int seed)
{
  /* generate second generation variants from seq of length seqlen.
     Use buffer varseq for variants.
     The original sequences came from seed */

  long matches = 0;

  /* make an exact copy */
  memcpy(varseq, seq, seqlen);

  /* substitutions */
  for(unsigned int i=0; i<seqlen; i++)
    {
      for (int v=1; v<5; v++)
        if (v != seq[i])
          {
            varseq[i] = v;
            if (hash_check_attach(varseq, seqlen, seed))
              matches++;
          }
      varseq[i] = seq[i];
    }

  /* deletions */
  if (seqlen > 1)
    memcpy(varseq, seq+1, seqlen-1);
  for(unsigned int i=0; i<seqlen; i++)
    {
      if ((i==0) || (seq[i] != seq[i-1]))
        {
          if (hash_check_attach(varseq, seqlen-1, seed))
            matches++;
        }
      varseq[i] = seq[i];
    }

  /* insertions */
  memcpy(varseq+1, seq, seqlen);
  for(unsigned int i=0; i<seqlen+1; i++)
    {
      for(int v=1; v<5; v++)
        {
          if((i==seqlen) || (v != seq[i]))
            {
              varseq[i] = v;
              if (hash_check_attach(varseq, seqlen+1, seed))
                matches++;
            }
        }
      if (i<seqlen)
        varseq[i] = seq[i];
    }
  return matches;
}

void fastidious_check_large_var(BloomFilter * bloom,
                                char * varseq,
                                char * buffer2,
                                int seed,
                                long * m,
                                long * v)
{
  /*
    bloom is a BloomFilter in which to enter the variants
    buffer1 is a buffer large enough to hold all sequences + 1 insertion
    buffer2 is a buffer large enough to hold all sequences + 2 insertions
    seed is the original seed
    m is where to store number of matches
    v is where to store number of variants
  */

  long variants = 0;
  long matches = 0;

  unsigned char * seq = (unsigned char*) db_getsequence(seed);
  unsigned long seqlen = db_getsequencelen(seed);

  /* make an exact copy */
  memcpy(varseq, seq, seqlen);

  /* substitutions */
  for(unsigned int i=0; i<seqlen; i++)
    {
      for (int v=1; v<5; v++)
        if (v != seq[i])
          {
            varseq[i] = v;
            variants++;
            if (bloom->get(varseq, seqlen))
              matches += fastidious_check_large_var_2(varseq,
                                                      seqlen,
                                                      buffer2,
                                                      seed);
          }
      varseq[i] = seq[i];
    }

  /* deletions */
  if (seqlen > 1)
    memcpy(varseq, seq+1, seqlen-1);
  for(unsigned int i=0; i<seqlen; i++)
    {
      if ((i==0) || (seq[i] != seq[i-1]))
        {
          variants++;
          if (bloom->get(varseq, seqlen-1))
            matches += fastidious_check_large_var_2(varseq,
                                                    seqlen-1,
                                                    buffer2,
                                                    seed);
        }
      varseq[i] = seq[i];
    }

  /* insertions */
  memcpy(varseq+1, seq, seqlen);
  for(unsigned int i=0; i<seqlen+1; i++)
    {
      for(int v=1; v<5; v++)
        {
          if((i==seqlen) || (v != seq[i]))
            {
              varseq[i] = v;
              variants++;
              if (bloom->get(varseq, seqlen+1))
                matches += fastidious_check_large_var_2(varseq,
                                                        seqlen+1,
                                                        buffer2,
                                                        seed);
            }
        }
      if (i<seqlen)
        varseq[i] = seq[i];
    }
  *m = matches;
  *v = variants;

#if 0
  long e = expected_variant_count((char*)seq, seqlen);
  if (variants != e)
    fprintf(logfile, "Incorrect number of variants: %ld Expected: %ld\n", variants, e);
#endif
}


int compare_amp(const void * a, const void * b)
{
  int * x = (int*) a;
  int * y = (int*) b;
  if (*x < *y)
    return -1;
  else if (*x > *y)
    return +1;
  else
    return 0;
}

static pthread_mutex_t light_mutex;
static long light_variants;
static long light_progress;
static long light_amplicon_count;
static int light_amplicon;
BloomFilter * bloomp;

void mark_light_thread(long t)
{
  (void) t;

  char * buffer1 = (char*) xmalloc(db_getlongestsequence() + 2);
  pthread_mutex_lock(&light_mutex);
  while (light_progress < light_amplicon_count)
    {
      int a = light_amplicon--;
      if (swarminfo[ampinfo[a].swarmid].mass < opt_boundary)
        {
          progress_update(++light_progress);
          pthread_mutex_unlock(&light_mutex);
          long v = fastidious_mark_small_var(bloomp, buffer1, a);
          pthread_mutex_lock(&light_mutex);
          light_variants += v;
        }
    }
  pthread_mutex_unlock(&light_mutex);
  free(buffer1);
}

static pthread_mutex_t heavy_mutex;
static long heavy_variants;
static long heavy_progress;
static long heavy_amplicon_count;
static int heavy_amplicon;
static long amplicons;

void check_heavy_thread(long t)
{
  (void) t;

  char * buffer1 = (char*) xmalloc(db_getlongestsequence() + 2);
  char * buffer2 = (char*) xmalloc(db_getlongestsequence() + 3);
  pthread_mutex_lock(&heavy_mutex);
  while ((heavy_amplicon < amplicons) && (heavy_progress < heavy_amplicon_count))
    {
      int a = heavy_amplicon++;
      if (swarminfo[ampinfo[a].swarmid].mass >= opt_boundary)
        {
          progress_update(++heavy_progress);
          pthread_mutex_unlock(&heavy_mutex);
          long m, v;
          fastidious_check_large_var(bloomp, buffer1, buffer2, a, &m, &v);
          pthread_mutex_lock(&heavy_mutex);
          heavy_variants += v;
        }
    }
  pthread_mutex_unlock(&heavy_mutex);
  free(buffer2);
  free(buffer1);
}

void algo_d1_run()
{
  unsigned long longestamplicon = db_getlongestsequence();
  amplicons = db_getsequencecount();

  threads_init();

  ampinfo = (struct ampinfo_s *)
    xmalloc(amplicons * sizeof(struct ampinfo_s));

  global_hits_alloc = longestamplicon * 7 + 4;
  global_hits_data = (int *) xmalloc(global_hits_alloc * sizeof(int));

  /* compute hash for all amplicons and store them in a hash table */

  hash_alloc(amplicons);

  progress_init("Hashing sequences:", amplicons);
  for(unsigned int i=0; i<amplicons; i++)
    {
      struct ampinfo_s * bp = ampinfo + i;
      bp->generation = 0;
      bp->swarmid = NO_SWARM;
      bp->next = NO_SWARM;
      bp->graft_cand = NO_SWARM;
      hash_insert(i);
      progress_update(i);
    }
  progress_done();

  unsigned char * dir = 0;
  unsigned long * hearray = 0;

  if (uclustfile)
    {
      dir = (unsigned char *) xmalloc(longestamplicon*longestamplicon);
      hearray = (unsigned long *)
        xmalloc(2 * longestamplicon * sizeof(unsigned long));
    }

  /* for each non-swarmed amplicon look for subseeds ... */
  long swarmid = 0;
  progress_init("Clustering:       ", amplicons);
  for(unsigned int seed = 0; seed < amplicons; seed++)
    {
      struct ampinfo_s * ap = ampinfo + seed;

      if (ap->swarmid == NO_SWARM)
        {
          /* start a new swarm with a new initial seed */

          ap->swarmid = swarmid;
          ap->generation = 0;
          ap->parent = NO_SWARM;
          ap->next = NO_SWARM;

          /* link up this initial seed in the list of swarms */
          current_swarm_tail = seed;

          /* initialize swarm stats */
          swarmsize = 0;
          swarm_maxgen = 0;
          abundance_sum = 0;
          singletons = 0;
          swarm_sumlen = 0;

          update_stats(seed);

          /* init list */
          global_hits_count = 0;

          /* find the first generation matches */
          process_seed(seed);

          /* sort hits */
          qsort(global_hits_data, global_hits_count,
                sizeof(int), compare_amp);

          /* add subseeds on list to current swarm */
          for(int i = 0; i < global_hits_count; i++)
            add_amp_to_swarm(global_hits_data[i]);

          /* find later generation matches */
          int subseed = ap->next;
          while(subseed != NO_SWARM)
            {
              /* process all subseeds of this generation */
              global_hits_count = 0;
              while(subseed != NO_SWARM)
                {
                  process_seed(subseed);
                  update_stats(subseed);
                  subseed = ampinfo[subseed].next;
                }

              /* sort all of this generation */
              qsort(global_hits_data, global_hits_count,
                    sizeof(int), compare_amp);

              /* add them to the swarm */
              for(int i = 0; i < global_hits_count; i++)
                add_amp_to_swarm(global_hits_data[i]);

              /* start with most abundant amplicon of next generation */
              if (global_hits_count)
                subseed = global_hits_data[0];
              else
                subseed = NO_SWARM;
            }

          if (swarmid >= swarminfo_alloc)
            {
              /* allocate memory for more swarms... */
              swarminfo_alloc += 1000;
              swarminfo =
                (struct swarminfo_s *) xrealloc (swarminfo,
                                                 swarminfo_alloc *
                                                 sizeof(swarminfo_s));
            }

          struct swarminfo_s * sp = swarminfo + swarmid;

          sp->seed = seed;
          sp->size = swarmsize;
          sp->mass = abundance_sum;
          sp->sumlen = swarm_sumlen;
          sp->singletons = singletons;
          sp->maxgen = swarm_maxgen;
          sp->last = current_swarm_tail;
          sp->attached = false;

          /* update overall stats */
          if (swarmsize > largest)
            largest = swarmsize;
          if (swarm_maxgen > maxgen)
            maxgen = swarm_maxgen;

          swarmid++;
        }
      progress_update(seed+1);
    }
  progress_done();

  long swarmcount = swarmid;

  swarmcount_adjusted = swarmcount;

  /* fastidious */

  if (opt_fastidious)
    {
      fprintf(logfile, "\n");
      fprintf(logfile, "Results before fastidious processing:\n");
      fprintf(logfile, "Number of swarms:  %ld\n", swarmcount);
      fprintf(logfile, "Largest swarm:     %d\n", largest);
      fprintf(logfile, "\n");

      long small_otus = 0;
      long amplicons_in_small_otus = 0;
      long nucleotides_in_small_otus = 0;

      progress_init("Counting amplicons in heavy and light swarms",
                    swarmcount);

      for(long i = 0; i < swarmcount; i++)
        {
          struct swarminfo_s * sp = swarminfo + i;
          if (sp->mass < opt_boundary)
            {
              amplicons_in_small_otus += sp->size;
              nucleotides_in_small_otus += sp->sumlen;
              small_otus++;
            }
          progress_update(i+1);
        }
      progress_done();

      long amplicons_in_large_otus = amplicons - amplicons_in_small_otus;
      long large_otus = swarmcount - small_otus;

      fprintf(logfile, "Heavy swarms: %ld, with %ld amplicons\n",
              large_otus, amplicons_in_large_otus);
      fprintf(logfile, "Light swarms: %ld, with %ld amplicons\n",
              small_otus, amplicons_in_small_otus);
      fprintf(logfile, "Total length of amplicons in light swarms: %ld\n",
              nucleotides_in_small_otus);

      if ((small_otus == 0) || (large_otus == 0))
        {
          fprintf(logfile, "Only light or heavy swarms found - "
                  "no need for further analysis.\n");
        }
      else
        {
          /* m: total size of Bloom filter in bits */
          /* k: number of hash functions */
          /* n: number of entries in the bloom filter */
          /* here: k=12 and m/n=18, that is 18 bits/entry */

          long bits = opt_bloom_bits; /* 18 */
          long k = int(bits * 0.693);    /* 12 */
          long m = bits * 7 * nucleotides_in_small_otus;

          long memtotal = arch_get_memtotal();
          long memused = arch_get_memused();

          if (opt_ceiling)
            {
              long memrest = 1024 * 1024 * opt_ceiling - memused;
              long new_bits = 8 * memrest / (7 * nucleotides_in_small_otus);
              if (new_bits < bits)
                {
                  if (new_bits < 2)
                    fatal("Insufficient memory remaining for Bloom filter");
                  fprintf(logfile, "Reducing memory used for Bloom filter due to --ceiling option.\n");
                  bits = new_bits;
                  k = int(bits * 0.693);
                  m = bits * 7 * nucleotides_in_small_otus;
                }
            }

          if (memused + m/8 > memtotal)
            {
              fprintf(logfile, "WARNING: Memory usage will probably exceed total amount of memory available.\n");
              fprintf(logfile, "Try to reduce memory footprint using the --bloom-bits or --ceiling options.\n");
            }

          fprintf(logfile,
                  "Bloom filter: bits=%ld, m=%ld, k=%ld, size=%.1fMB\n",
                  bits, m, k, 1.0 * m / (8*1024*1024));

          bloomp = new BloomFilter(m, k);
          char * buffer1 = (char*) xmalloc(db_getlongestsequence() + 2);
          char * buffer2 = (char*) xmalloc(db_getlongestsequence() + 3);

          progress_init("Adding light swarm amplicons to Bloom filter",
                        amplicons_in_small_otus);

          /* process amplicons in order from least to most abundant */
          /* but stop when all amplicons in small otus are processed */

          light_variants = 0;

#if 1
          pthread_mutex_init(&light_mutex, NULL);
          light_progress = 0;
          light_amplicon_count = amplicons_in_small_otus;
          light_amplicon = amplicons - 1;
          ThreadRunner * tr = new ThreadRunner(opt_threads, mark_light_thread);
          tr->run();
          delete tr;
          pthread_mutex_destroy(&light_mutex);
#else
          int a = amplicons - 1;
          long x = 0;
          while (x < amplicons_in_small_otus)
            {
              if (swarminfo[ampinfo[a].swarmid].mass < opt_boundary)
                {
                  light_variants += fastidious_mark_small_var(bloomp, buffer1, a);
                  x++;
                  progress_update(x);
                }
              a--;
            }
#endif

          progress_done();

          fprintf(logfile,
                  "Generated %ld variants from light swarms\n", light_variants);

          progress_init("Checking heavy swarm amplicons against Bloom filter",
                        amplicons_in_large_otus);

          /* process amplicons in order from most to least abundant */
          /* but stop when all amplicons in large otus are processed */

          pthread_mutex_init(&graft_mutex, NULL);

          heavy_variants = 0;

#if 1
          pthread_mutex_init(&heavy_mutex, NULL);
          heavy_progress = 0;
          heavy_amplicon_count = amplicons_in_large_otus;
          heavy_amplicon = 0;
          ThreadRunner * heavy_tr = new ThreadRunner(opt_threads,
                                                     check_heavy_thread);
          heavy_tr->run();
          delete heavy_tr;
          pthread_mutex_destroy(&heavy_mutex);
#else
          long i = 0;

          for(int a = 0; (a < amplicons) && (i < amplicons_in_large_otus); a++)
            {
              int swarmid = ampinfo[a].swarmid;
              int mass = swarminfo[swarmid].mass;
              if (mass >= opt_boundary)
                {
                  long m, v;
                  fastidious_check_large_var(bloomp, buffer1, buffer2,
                                             a, &m, &v);
                  heavy_variants += v;
                  progress_update(++i);
                }
            }
#endif

          progress_done();

          free(buffer1);
          free(buffer2);

          delete bloomp;

          pthread_mutex_destroy(&graft_mutex);

          fprintf(logfile, "Heavy variants: %ld\n", heavy_variants);
          fprintf(logfile, "Got %ld graft candidates\n", graft_candidates);
          int grafts = attach_candidates(amplicons);
          fprintf(logfile, "Made %d grafts\n", grafts);
          fprintf(logfile, "\n");
        }
    }


  /* dump swarms */

  progress_init("Writing swarms:   ", swarmcount);

  if (opt_mothur)
    fprintf(outfile, "swarm_%ld\t%ld", opt_differences, swarmcount_adjusted);

  for(int i = 0; i < swarmcount; i++)
    {
      if (!swarminfo[i].attached)
        {
          int seed = swarminfo[i].seed;
          for (int a = seed;
               a >= 0;
               a = ampinfo[a].next)
            {
              if (opt_mothur)
                {
                  if (a == seed)
                    fputc('\t', outfile);
                  else
                    fputc(',', outfile);
                }
              else
                {
                  if (a != seed)
                    fputc(SEPCHAR, outfile);
                }
              fprint_id(outfile, a);
            }
          if (!opt_mothur)
            fputc('\n', outfile);
        }
      progress_update(i+1);
    }

  if (opt_mothur)
    fputc('\n', outfile);

  progress_done();


  /* dump seeds in fasta format with sum of abundances */

  if (opt_seeds)
    {
      progress_init("Writing seeds:    ", swarmcount);
      for(int i=0; i < swarmcount; i++)
        {
          if (!swarminfo[i].attached)
            {
              int seed = swarminfo[i].seed;
              fprintf(fp_seeds, ">");
              fprint_id_with_new_abundance(fp_seeds, seed, swarminfo[i].mass);
              fprintf(fp_seeds, "\n");
              db_fprintseq(fp_seeds, seed, 0);
            }
          progress_update(i+1);
        }
      progress_done();
    }


  /* output internal structure */

  if (opt_internal_structure)
    {
      unsigned int cluster_no = 0;

      progress_init("Writing structure:", swarmcount);

      for(unsigned int swarmid = 0; swarmid < swarmcount ; swarmid++)
        {
          if (!swarminfo[swarmid].attached)
            {
              int seed = swarminfo[swarmid].seed;

              struct ampinfo_s * bp = ampinfo + seed;

              for (int a = bp->next;
                   a >= 0;
                   a = ampinfo[a].next)
                {
                  long graft_parent = ampinfo[a].graft_cand;
                  if (graft_parent != NO_SWARM)
                    {
                      fprint_id_noabundance(internal_structure_file, graft_parent);
                      fprintf(internal_structure_file, "\t");
                      fprint_id_noabundance(internal_structure_file, a);
                      fprintf(internal_structure_file,
                              "\t%d\t%d\t%d\n",
                              2,
                              cluster_no + 1,
                              ampinfo[graft_parent].generation + 1);
                    }

                  long parent = ampinfo[a].parent;
                  if (parent != NO_SWARM)
                    {
                      int diff = 1;
                      if (duplicates_found)
                        {
                          unsigned long parentseqlen = db_getsequencelen(parent);
                          unsigned long ampseqlen = db_getsequencelen(a);
                          if (parentseqlen == ampseqlen)
                            {
                              unsigned char * parentseq = (unsigned char *) db_getsequence(parent);
                              unsigned char * ampseq = (unsigned char *) db_getsequence(a);
                              if (memcmp(parentseq, ampseq, parentseqlen) == 0)
                                diff = 0;
                            }
                        }
                      fprint_id_noabundance(internal_structure_file, parent);
                      fprintf(internal_structure_file, "\t");
                      fprint_id_noabundance(internal_structure_file, a);
                      fprintf(internal_structure_file,
                              "\t%d\t%d\t%d\n",
                              diff,
                              cluster_no + 1,
                              ampinfo[a].generation);
                    }
                }

              cluster_no++;
            }
          progress_update(swarmid);
        }
      progress_done();
    }


  /* output swarm in uclust format */

  if (uclustfile)
    {
      unsigned int cluster_no = 0;

      progress_init("Writing UCLUST:   ", swarmcount);

      for(unsigned int swarmid = 0; swarmid < swarmcount ; swarmid++)
        {
          if (!swarminfo[swarmid].attached)
            {
              int seed = swarminfo[swarmid].seed;

              struct ampinfo_s * bp = ampinfo + seed;

              fprintf(uclustfile, "C\t%u\t%d\t*\t*\t*\t*\t*\t",
                      cluster_no,
                      swarminfo[swarmid].size);
              fprint_id(uclustfile, seed);
              fprintf(uclustfile, "\t*\n");

              fprintf(uclustfile, "S\t%u\t%lu\t*\t*\t*\t*\t*\t",
                      cluster_no,
                      db_getsequencelen(seed));
              fprint_id(uclustfile, seed);
              fprintf(uclustfile, "\t*\n");

              for (int a = bp->next;
                   a >= 0;
                   a = ampinfo[a].next)
                {
                  char * dseq = db_getsequence(a);
                  char * dend = dseq + db_getsequencelen(a);
                  char * qseq = db_getsequence(seed);
                  char * qend = qseq + db_getsequencelen(seed);

                  unsigned long nwscore = 0;
                  unsigned long nwdiff = 0;
                  char * nwalignment = NULL;
                  unsigned long nwalignmentlength = 0;

                  nw(dseq, dend, qseq, qend,
                     score_matrix_63, penalty_gapopen, penalty_gapextend,
                     & nwscore, & nwdiff, & nwalignmentlength, & nwalignment,
                     dir, hearray, 0, 0);

                  double percentid = 100.0 * (nwalignmentlength - nwdiff) /
                    nwalignmentlength;

                  fprintf(uclustfile,
                          "H\t%d\t%lu\t%.1f\t+\t0\t0\t%s\t",
                          cluster_no,
                          db_getsequencelen(a),
                          percentid,
                          nwdiff > 0 ? nwalignment : "=");

                  fprint_id(uclustfile, a);
                  fprintf(uclustfile, "\t");
                  fprint_id(uclustfile, seed);
                  fprintf(uclustfile, "\n");

                  if (nwalignment)
                    free(nwalignment);
                }

              cluster_no++;
            }
          progress_update(swarmid);
        }
      progress_done();
    }

  /* output statistics to file */

  if (statsfile)
    {
      progress_init("Writing stats:    ", swarmcount);
      for(long i = 0; i < swarmcount; i++)
        {
          swarminfo_s * sp = swarminfo + i;
          if (!sp->attached)
            {
              fprintf(statsfile, "%d\t%ld\t", sp->size, sp->mass);
              fprint_id_noabundance(statsfile, sp->seed);
              fprintf(statsfile, "\t%lu\t%d\t%d\t%d\n",
                      db_getabundance(sp->seed),
                      sp->singletons, sp->maxgen, sp->maxgen);
            }
          progress_update(i);
        }
      progress_done();
    }


  fprintf(logfile, "\n");
  fprintf(logfile, "Number of swarms:  %ld\n", swarmcount_adjusted);
  fprintf(logfile, "Largest swarm:     %d\n", largest);
  fprintf(logfile, "Max generations:   %d\n", maxgen);

  threads_done();

  hash_free();

  if(swarminfo)
    free(swarminfo);

  free(ampinfo);

  free(global_hits_data);

  if (uclustfile)
    {
      free(dir);
      free(hearray);
    }

#ifdef HASHSTATS
  fprintf(logfile, "Tries:      %12lu\n", tries);
  fprintf(logfile, "Bloom m:    %12lu\n", bloom_matches);
  fprintf(logfile, "Hits:       %12lu\n", hits);
  fprintf(logfile, "Success:    %12lu\n", success);
  fprintf(logfile, "Bingo:      %12lu\n", bingo);
  fprintf(logfile, "Collisions: %12lu\n", collisions);
#endif
}
