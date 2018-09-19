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

/* Information about each amplicon */

static struct ampinfo_s
{
  int swarmid;
  int parent;
  int generation;
  int next;       /* amp id of next amplicon in swarm */
  int graft_cand; /* amp id of potential grafting parent (fastidious) */
  int link_start;
  int link_count;
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

static struct graft_cand
{
  int parent;
  int child;
} * graft_array = 0;

/* Information about potential grafts */
static long graft_candidates = 0;
static pthread_mutex_t graft_mutex;

#define NO_SWARM (-1)

static int current_swarm_tail = 0;

static long swarminfo_alloc = 0;

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

static int * global_hits_data = 0;
static int global_hits_alloc = 0;
static int global_hits_count = 0;

static unsigned long longestamplicon = 0;

static long amplicons = 0;

static pthread_mutex_t heavy_mutex;
static long heavy_variants = 0;
static long heavy_progress = 0;
static long heavy_amplicon_count = 0;
static int heavy_amplicon = 0;

static pthread_mutex_t light_mutex;
static long light_variants = 0;
static long light_progress = 0;
static long light_amplicon_count = 0;
static int light_amplicon = 0;

static unsigned long network_alloc = 1024 * 1024;
static int * network = 0;
static unsigned long network_count = 0;
static pthread_mutex_t network_mutex;
static long network_amp = 0;

static struct bloom_s * bloom_a = 0; // Bloom filter for amplicons

struct bloomflex_s * bloom_f = 0; // Huge Bloom filter for fastidious

inline bool check_amp_identical(unsigned int amp1,
                                unsigned int amp2)
{
  unsigned int amp1_seqlen = db_getsequencelen(amp1);
  unsigned int amp2_seqlen = db_getsequencelen(amp2);
  
  if (amp1_seqlen == amp2_seqlen)
    {
      char * amp1_sequence = db_getsequence(amp1);
      char * amp2_sequence = db_getsequence(amp2);
      
      return memcmp(amp1_sequence,
                    amp2_sequence,
                    nt_bytelength(amp1_seqlen)) == 0;
    }
  else
    return false;
}

inline void hash_insert(int amp)
{
  /* find the first empty bucket */
  unsigned long hash = db_gethash(amp);
  unsigned int j = hash_getindex(hash);
  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, hash) &&
          check_amp_identical(amp, hash_get_data(j)))
        duplicates_found++;
      j = hash_getnextindex(j);
    }
  
  hash_set_occupied(j);
  hash_set_value(j, hash);
  hash_set_data(j, amp);

  bloom_set(bloom_a, hash);
}


/******************** FASTIDIOUS START ********************/


void attach(int seed, int amp)
{
  /* graft light swarm (amp) on heavy swarm (seed) */

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
                       unsigned int seqlen,
                       struct var_s * var,
                       int seed)
{
  /* seed is the original large swarm seed */

  /* compute hash and corresponding hash table index */
  unsigned long hash = var->hash;
  unsigned int j = hash_getindex(hash);

  /* find matching buckets */

  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, hash))
        {
          /* check that mass is below threshold */
          int amp = hash_get_data(j);

          /* make absolutely sure sequences are identical */
          char * ampseq = db_getsequence(amp);
          unsigned long ampseqlen = db_getsequencelen(amp);
          if (check_variant(seq, seqlen, var, ampseq, ampseqlen))
            {
              add_graft_candidate(seed, amp);
              return true;
            }
        }
      j = hash_getnextindex(j);
    }
  return false;
}

inline long check_heavy_var_2(char * seq,
                              size_t seqlen,
                              int seed,
                              struct var_s * variant_list)
{
  /* Check second generation microvariants of the heavy swarm amplicons
     and see if any of them are identical to a light swarm amplicon. */

  long matches = 0;
  unsigned int variant_count = 0;

  unsigned long hash = zobrist_hash((unsigned char*)seq, seqlen);
  generate_variants(seq, seqlen, hash,
                    variant_list, & variant_count, false);

  for(unsigned int i=0; i < variant_count; i++)
    if (bloom_get(bloom_a, variant_list[i].hash) &&
        hash_check_attach(seq, seqlen, variant_list + i, seed))
      matches++;

  return matches;
}

void check_heavy_var(struct bloomflex_s * bloom,
                     char * varseq,
                     int seed,
                     long * m,
                     long * v,
                     struct var_s * variant_list,
                     struct var_s * variant_list2)
{
  /*
    bloom is a bloom filter in which to check the variants
    varseq is a buffer large enough to hold all sequences + 1 insertion
    seed is the original seed
    m is where to store number of matches
    v is where to store number of variants
    variant_list and variant_list2 are lists to hold the 1st and 2nd
    generation of microvariants
  */

  /*
    Generate microvariants of the heavy swarm amplicons, forming
    "virtual" amplicons. Check with the bloom filter if any
    of these are identical to the microvariants of the
    light swarm amplicons. If there is a match we have a potential
    link. To find which light amplicon it could link to, we have
    to generate the second generation microvariants and check
    these against the light swarm amplicons.
  */

  unsigned int variant_count = 0;
  long matches = 0;

  char * sequence = db_getsequence(seed);
  unsigned int seqlen = db_getsequencelen(seed);
  unsigned long hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash,
                    variant_list, & variant_count, false);

  for(unsigned int i = 0; i < variant_count; i++)
    {
      struct var_s * var = variant_list + i;
      if (bloomflex_get(bloom, var->hash))
        {
          int varlen = 0;
          generate_variant_sequence(sequence, seqlen,
                                    var, varseq, & varlen);
          matches += check_heavy_var_2(varseq,
                                       varlen,
                                       seed,
                                       variant_list2);
        }
    }

  *m = matches;
  *v = variant_count;
}

void check_heavy_thread(long t)
{
  (void) t;

  struct var_s * variant_list
    = (struct var_s *) xmalloc(sizeof(struct var_s) *
                               (7 * longestamplicon + 4));
  struct var_s * variant_list2
    = (struct var_s *) xmalloc(sizeof(struct var_s) *
                               (7 * (longestamplicon+1) + 4));

  char * buffer1 = (char*) xmalloc(db_getlongestsequence() + 2);
  pthread_mutex_lock(&heavy_mutex);
  while ((heavy_amplicon < amplicons) &&
         (heavy_progress < heavy_amplicon_count))
    {
      int a = heavy_amplicon++;
      if (swarminfo[ampinfo[a].swarmid].mass >= opt_boundary)
        {
          progress_update(++heavy_progress);
          pthread_mutex_unlock(&heavy_mutex);
          long m, v;
          check_heavy_var(bloom_f, buffer1, a, &m, &v,
                          variant_list, variant_list2);
          pthread_mutex_lock(&heavy_mutex);
          heavy_variants += v;
        }
    }
  pthread_mutex_unlock(&heavy_mutex);
  free(buffer1);

  free(variant_list2);
  free(variant_list);
}

long mark_light_var(struct bloomflex_s * bloom,
                    int seed,
                    struct var_s * variant_list)
{
  /*
    add all microvariants of seed to Bloom filter

    bloom is a BloomFilter in which to enter the variants
    seed is the original seed
  */

  hash_insert(seed);

  unsigned int variant_count = 0;

  char * sequence = db_getsequence(seed);
  unsigned int seqlen = db_getsequencelen(seed);
  unsigned long hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash,
                    variant_list, & variant_count, false);

  for(unsigned int i = 0; i < variant_count; i++)
    bloomflex_set(bloom, variant_list[i].hash);

  return variant_count;
}

void mark_light_thread(long t)
{
  (void) t;

  struct var_s * variant_list
    = (struct var_s *) xmalloc(sizeof(struct var_s) *
                               (7 * longestamplicon + 4));

  pthread_mutex_lock(&light_mutex);
  while (light_progress < light_amplicon_count)
    {
      int a = light_amplicon--;
      if (swarminfo[ampinfo[a].swarmid].mass < opt_boundary)
        {
          progress_update(++light_progress);
          pthread_mutex_unlock(&light_mutex);
          long v = mark_light_var(bloom_f, a, variant_list);
          pthread_mutex_lock(&light_mutex);
          light_variants += v;
        }
    }
  pthread_mutex_unlock(&light_mutex);

  free(variant_list);
}


/******************** FASTIDIOUS END ********************/


inline void find_variant_matches(int seed,
                                 var_s * var,
                                 int * hits_data,
                                 int * hits_count)
{
  if (! bloom_get(bloom_a, var->hash))
    return;

  /* compute hash and corresponding hash table index */

  unsigned int j = hash_getindex(var->hash);

  /* find matching buckets */

  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, var->hash))
        {
          int amp = hash_get_data(j);

          /* avoid self */
          if (seed != amp)
            if (opt_no_otu_breaking ||
                (db_getabundance(seed) >= db_getabundance(amp)))
              {
                char * seed_sequence = db_getsequence(seed);
                unsigned int seed_seqlen = db_getsequencelen(seed);

                char * amp_sequence = db_getsequence(amp);
                unsigned int amp_seqlen = db_getsequencelen(amp);

                if (check_variant(seed_sequence, seed_seqlen,
                                  var,
                                  amp_sequence, amp_seqlen))
                  {
                    hits_data[(*hits_count)++] = amp;
                    break;
                  }
              }
        }
      j = hash_getnextindex(j);
    }
}

void check_variants(int seed,
                    var_s * variant_list,
                    int * hits_data,
                    int * hits_count)
{
  unsigned int variant_count = 0;
  * hits_count = 0;

  char * sequence = db_getsequence(seed);
  unsigned int seqlen = db_getsequencelen(seed);
  unsigned long hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash,
                    variant_list, & variant_count, true);

  for(unsigned int i = 0; i < variant_count; i++)
    find_variant_matches(seed, variant_list + i, hits_data, hits_count);
}

void network_thread(long t)
{
  (void) t;

  int hits_count = 0;
  int * hits_data
    = (int *) xmalloc((7 * longestamplicon + 5) * sizeof(int));

  struct var_s * variant_list
    = (var_s *) xmalloc((7 * longestamplicon + 5) * sizeof(struct var_s));

  pthread_mutex_lock(&network_mutex);
  while (network_amp < amplicons)
    {
      int amp = network_amp++;
      progress_update(amp);

      pthread_mutex_unlock(&network_mutex);

      hits_count = 0;
      check_variants(amp, variant_list, hits_data, & hits_count);
      pthread_mutex_lock(&network_mutex);

      ampinfo[amp].link_start = network_count;
      ampinfo[amp].link_count = hits_count;

      if (network_count + hits_count > network_alloc)
        {
          while (network_count + hits_count > network_alloc)
            network_alloc += 1024 * 1024;
          network = (int*) xrealloc(network, network_alloc * sizeof(int));
        }

      for(int i=0; i < hits_count; i++)
        network[network_count++] = hits_data[i];
    }
  pthread_mutex_unlock(&network_mutex);

  free(variant_list);
  free(hits_data);
}

void process_seed(int seed)
{
  /* update swarm stats */
  struct ampinfo_s * bp = ampinfo + seed;

  swarmsize++;
  if (bp->generation > swarm_maxgen)
    swarm_maxgen = bp->generation;
  unsigned long abundance = db_getabundance(seed);
  abundance_sum += abundance;
  if (abundance == 1)
    singletons++;
  swarm_sumlen += db_getsequencelen(seed);

  int s = ampinfo[seed].link_start;
  int c = ampinfo[seed].link_count;

  if (global_hits_count + c > global_hits_alloc)
    {
      while (global_hits_count + c > global_hits_alloc)
        global_hits_alloc += 4096;
      global_hits_data = (int*)xrealloc(global_hits_data,
                                        global_hits_alloc * sizeof(int));
    }

  for(int i = 0; i < c; i++)
    {
      int amp = network[s + i];

      if (ampinfo[amp].swarmid == NO_SWARM)
        {
          global_hits_data[global_hits_count++] = amp;

          /* update info */
          ampinfo[amp].swarmid = ampinfo[seed].swarmid;
          ampinfo[amp].generation = ampinfo[seed].generation + 1;
          ampinfo[amp].parent = seed;
        }
    }
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

inline void add_amp_to_swarm(int amp)
{
  /* add to swarm */
  ampinfo[current_swarm_tail].next = amp;
  current_swarm_tail = amp;
}

void algo_d1_run()
{
  longestamplicon = db_getlongestsequence();
  amplicons = db_getsequencecount();

  ampinfo = (struct ampinfo_s *)
    xmalloc(amplicons * sizeof(struct ampinfo_s));

  global_hits_alloc = longestamplicon * 7 + 4 + 1;
  global_hits_data = (int *) xmalloc(global_hits_alloc * sizeof(int));

  /* compute hash for all amplicons and store them in a hash table */

  hash_alloc(amplicons);
  bloom_a = bloom_init(hash_get_tablesize());

  duplicates_found = 0;

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

  if (duplicates_found)
    {
      fprintf(logfile,
              "WARNING: %lu duplicated sequences detected.\n"
              "Please consider dereplicating your data for optimal results.\n",
              duplicates_found);
    }

  unsigned char * dir = 0;
  unsigned long * hearray = 0;

  if (uclustfile)
    {
      dir = (unsigned char *) xmalloc(longestamplicon*longestamplicon);
      hearray = (unsigned long *)
        xmalloc(2 * longestamplicon * sizeof(unsigned long));
    }

  /* for all amplicons, generate list of matching amplicons */

  network = (int*) xmalloc(network_alloc * sizeof(int));
  network_count = 0;

  pthread_mutex_init(&network_mutex, NULL);
  network_amp = 0;
  progress_init("Building network: ", amplicons);
  ThreadRunner * network_tr = new ThreadRunner(opt_threads, network_thread);
  network_tr->run();
  delete network_tr;
  pthread_mutex_destroy(&network_mutex);

  progress_done();

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
              swarminfo_alloc += 1024;
              swarminfo =
                (struct swarminfo_s *) xrealloc(swarminfo,
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

  free(global_hits_data);

  free(network);
  network = 0;

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
          /* here: k=11 and m/n=18, that is 16 bits/entry */

          long bits = opt_bloom_bits; /* 16 */
          // long k = int(bits * 0.693);    /* 11 */
          long k = int(bits * 0.4);    /* 6 */
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
                  // k = int(bits * 0.693);
                  k = int(bits * 0.4);
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


          bloom_f = bloomflex_init(m/8, k);


          /* Empty the old hash and bloom filter
             before we reinsert only the light swarm amplicons */

          hash_zap();
          bloom_zap(bloom_a);

          progress_init("Adding light swarm amplicons to Bloom filter",
                        amplicons_in_small_otus);

          /* process amplicons in order from least to most abundant */
          /* but stop when all amplicons in small otus are processed */

          light_variants = 0;

          pthread_mutex_init(&light_mutex, NULL);
          light_progress = 0;
          light_amplicon_count = amplicons_in_small_otus;
          light_amplicon = amplicons - 1;
          ThreadRunner * tr = new ThreadRunner(opt_threads, mark_light_thread);
          tr->run();
          delete tr;
          pthread_mutex_destroy(&light_mutex);

          progress_done();

          fprintf(logfile,
                  "Generated %ld variants from light swarms\n",
                  light_variants);

          progress_init("Checking heavy swarm amplicons against Bloom filter",
                        amplicons_in_large_otus);

          /* process amplicons in order from most to least abundant */
          /* but stop when all amplicons in large otus are processed */

          pthread_mutex_init(&graft_mutex, NULL);

          heavy_variants = 0;

          pthread_mutex_init(&heavy_mutex, NULL);
          heavy_progress = 0;
          heavy_amplicon_count = amplicons_in_large_otus;
          heavy_amplicon = 0;
          ThreadRunner * heavy_tr = new ThreadRunner(opt_threads,
                                                     check_heavy_thread);
          heavy_tr->run();
          delete heavy_tr;
          pthread_mutex_destroy(&heavy_mutex);

          progress_done();

          bloomflex_exit(bloom_f);

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
                      fprint_id_noabundance(internal_structure_file,
                                            graft_parent);
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

  bloom_exit(bloom_a);
  hash_free();

  if (swarminfo)
    free(swarminfo);

  free(ampinfo);

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
