/*
    SWARM

    Copyright (C) 2012-2020 Torbjorn Rognes and Frederic Mahe

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

constexpr unsigned int one_kilobyte {1 << 10};
constexpr unsigned int one_megabyte {one_kilobyte * one_kilobyte};

/* Information about each amplicon */

static struct ampinfo_s
{
  unsigned int swarmid;
  unsigned int parent;
  unsigned int generation;
  unsigned int next;       /* amp id of next amplicon in swarm */
  unsigned int graft_cand; /* amp id of potential grafting parent (fastid.) */
  unsigned int link_start;
  unsigned int link_count;
} * ampinfo = nullptr;

/* Information about each swarm (OTU) */

static struct swarminfo_s
{
  uint64_t mass; /* the sum of abundances of amplicons in this swarm */
  uint64_t sumlen; /* sum of length of amplicons in swarm */
  unsigned int seed; /* amplicon id of the initial seed of this swarm */
  unsigned int last; /* amplicon id of the last seed in this swarm */
  unsigned int size; /* total number of amplicons in this swarm */
  unsigned int singletons; /* number of amplicons with abundance 1 */
  unsigned int maxgen; /* the generation of the amplicon farthest from seed */
  bool attached; /* this is a small swarm attached to a large (fastidious) */
  char dummy[3]; /* alignment padding only */
} * swarminfo = nullptr;

static struct graft_cand
{
  unsigned int parent;
  unsigned int child;
} * graft_array = nullptr;

/* Information about potential grafts */
static int64_t graft_candidates = 0;
static pthread_mutex_t graft_mutex;

constexpr unsigned int no_swarm {UINT_MAX};

static unsigned int current_swarm_tail = 0;

static uint64_t swarminfo_alloc = 0;

/* overall statistics */
static unsigned int maxgen = 0;
static unsigned int largest = 0;
static uint64_t swarmcount_adjusted = 0;

/* per swarm statistics */
static unsigned int singletons = 0;
static uint64_t abundance_sum = 0; /* = mass */
static unsigned int swarmsize = 0;
static unsigned int swarm_maxgen = 0;
static uint64_t swarm_sumlen = 0;

static unsigned int * global_hits_data = nullptr;
static unsigned int global_hits_alloc = 0;
static unsigned int global_hits_count = 0;

static unsigned int longestamplicon = 0;

static unsigned int amplicons = 0;

static pthread_mutex_t heavy_mutex;
static uint64_t heavy_variants = 0;
static uint64_t heavy_progress = 0;
static uint64_t heavy_amplicon_count = 0;
static unsigned int heavy_amplicon = 0;

static pthread_mutex_t light_mutex;
static uint64_t light_variants = 0;
static uint64_t light_progress = 0;
static uint64_t light_amplicon_count = 0;
static unsigned int light_amplicon = 0;

static uint64_t network_alloc = one_megabyte;
static unsigned int * network = nullptr;
static unsigned int network_count = 0;
static pthread_mutex_t network_mutex;
static unsigned int network_amp = 0;

static struct bloom_s * bloom_a = nullptr; // Bloom filter for amplicons

static struct bloomflex_s * bloom_f = nullptr; // Huge Bloom filter for fastidious

void attach(unsigned int seed, unsigned int amp);
void add_graft_candidate(unsigned int seed, unsigned int amp);
int compare_grafts(const void * a, const void * b);
unsigned int attach_candidates(unsigned int amplicon_count);
bool hash_check_attach(char * seq,
                       unsigned int seqlen,
                       struct var_s * var,
                       unsigned int seed);
void check_heavy_var(struct bloomflex_s * bloom,
                     char * varseq,
                     unsigned int seed,
                     uint64_t * m,
                     uint64_t * v,
                     struct var_s * variant_list,
                     struct var_s * variant_list2);
void check_heavy_thread(int64_t t);
uint64_t mark_light_var(struct bloomflex_s * bloom,
                        unsigned int seed,
                        struct var_s * variant_list);
void mark_light_thread(int64_t t);
void check_variants(unsigned int seed,
                    var_s * variant_list,
                    unsigned int * hits_data,
                    unsigned int * hits_count);
void network_thread(int64_t t);
void process_seed(unsigned int seed);
int compare_amp(const void * a, const void * b);
int compare_mass(const void * a, const void * b);


inline bool check_amp_identical(unsigned int amp1,
                                unsigned int amp2)
{
  /* amplicon are identical if they have the same length, and the
     exact same sequence */
  unsigned int amp1_seqlen = db_getsequencelen(amp1);
  unsigned int amp2_seqlen = db_getsequencelen(amp2);
  auto status {false};

  if (amp1_seqlen == amp2_seqlen) {
    char * amp1_sequence = db_getsequence(amp1);
    char * amp2_sequence = db_getsequence(amp2);

    if (memcmp(amp1_sequence,
               amp2_sequence,
               nt_bytelength(amp1_seqlen)) == 0) {
      status = true;
    }
  }

  return status;
}

inline void hash_insert(unsigned int amp)
{
  /* find the first empty bucket */
  uint64_t hash = db_gethash(amp);
  uint64_t j = hash_getindex(hash);
  bool duplicate = false;
  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, hash) &&
          check_amp_identical(amp, hash_get_data(j)))
        duplicate = true;
      j = hash_getnextindex(j);
    }

  if (duplicate)
    duplicates_found++;

  hash_set_occupied(j);
  hash_set_value(j, hash);
  hash_set_data(j, amp);

  bloom_set(bloom_a, hash);
}


/******************** FASTIDIOUS START ********************/


void attach(unsigned int seed, unsigned int amp)
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

void add_graft_candidate(unsigned int seed, unsigned int amp)
{
  pthread_mutex_lock(&graft_mutex);
  graft_candidates++;
  if ((ampinfo[amp].graft_cand == no_swarm)||(ampinfo[amp].graft_cand > seed))
    ampinfo[amp].graft_cand = seed;
  pthread_mutex_unlock(&graft_mutex);
}

int compare_grafts(const void * a, const void * b)
{
  const auto * x = static_cast<const struct graft_cand *>(a);
  const auto * y = static_cast<const struct graft_cand *>(b);
  auto status {0};

  /* replace with a three-way comparison '<=>' in a few years */
  // status = x->parent <=> y->parent : -1 : 0 : +1;
  // if (status == 0) {
  //   status = x->child <=> y->child : -1 : 0 : +1;
  // }

  if (x->parent < y->parent) {
    status = -1;
  }
  else if (x->parent > y->parent) {
    status = +1;
  }
  else {
    if (x->child < y->child) {
      status = -1;
    }
    else if (x->child > y->child) {
      status = +1;
    }
    else {
      status = 0;
    }
  }

  return status;
}

unsigned int attach_candidates(unsigned int amplicon_count)
{
  /* count pairs */
  unsigned int pair_count = 0;
  for(auto i = 0U; i < amplicon_count; i++)
    if (ampinfo[i].graft_cand != no_swarm)
      pair_count++;

  unsigned int grafts = 0;
  progress_init("Grafting light swarms on heavy swarms", pair_count);

  /* allocate memory */
  graft_array = static_cast<struct graft_cand *>
    (xmalloc(pair_count * sizeof(struct graft_cand)));

  /* fill in */
  unsigned int j = 0;
  for(auto i = 0U; i < amplicon_count; i++)
    if (ampinfo[i].graft_cand != no_swarm)
      {
        graft_array[j].parent = ampinfo[i].graft_cand;
        graft_array[j].child = i;
        j++;
      }

  /* sort */
  qsort(graft_array, pair_count, sizeof(struct graft_cand), compare_grafts);

  /* attach in order */
  for(auto i = 0U; i < pair_count; i++)
    {
      unsigned int parent = graft_array[i].parent;
      unsigned int child  = graft_array[i].child;

      if (swarminfo[ampinfo[child].swarmid].attached)
        {
          /* this light swarm is already attached */
          ampinfo[child].graft_cand = no_swarm;
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
  xfree(graft_array);
  return grafts;
}

bool hash_check_attach(char * seq,
                       unsigned int seqlen,
                       struct var_s * var,
                       unsigned int seed)
{
  /* seed is the original large swarm seed */

  /* compute hash and corresponding hash table index */
  uint64_t hash = var->hash;
  uint64_t j = hash_getindex(hash);

  /* find matching buckets */

  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, hash))
        {
          /* check that mass is below threshold */
          unsigned int amp = hash_get_data(j);

          /* make absolutely sure sequences are identical */
          char * ampseq = db_getsequence(amp);
          unsigned int ampseqlen = db_getsequencelen(amp);
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

inline uint64_t check_heavy_var_2(char * seq,
                                  unsigned int seqlen,
                                  unsigned int seed,
                                  struct var_s * variant_list)
{
  /* Check second generation microvariants of the heavy swarm amplicons
     and see if any of them are identical to a light swarm amplicon. */

  uint64_t matches = 0;
  unsigned int variant_count = 0;

  uint64_t hash = zobrist_hash(reinterpret_cast<unsigned char *>(seq), seqlen);
  generate_variants(seq, seqlen, hash, variant_list, & variant_count);

  for(auto i = 0U; i < variant_count; i++)
    if (bloom_get(bloom_a, variant_list[i].hash) &&
        hash_check_attach(seq, seqlen, variant_list + i, seed))
      matches++;

  return matches;
}

void check_heavy_var(struct bloomflex_s * bloom,
                     char * varseq,
                     unsigned int seed,
                     uint64_t * m,
                     uint64_t * v,
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
  uint64_t matches = 0;

  char * sequence = db_getsequence(seed);
  unsigned int seqlen = db_getsequencelen(seed);
  uint64_t hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash, variant_list, & variant_count);

  for(auto i = 0U; i < variant_count; i++)
    {
      struct var_s * var = variant_list + i;
      if (bloomflex_get(bloom, var->hash))
        {
          unsigned int varlen = 0;
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

void check_heavy_thread(int64_t t)
{
  (void) t;

  struct var_s * variant_list = static_cast<struct var_s *>
    (xmalloc(sizeof(struct var_s) * (7 * longestamplicon + 4)));
  struct var_s * variant_list2 = static_cast<struct var_s *>
    (xmalloc(sizeof(struct var_s) * (7 * (longestamplicon+1) + 4)));

  size_t size = 8 * ((db_getlongestsequence() + 2 + 31) / 32);
  char * buffer1 = static_cast<char *>(xmalloc(size));
  pthread_mutex_lock(&heavy_mutex);
  while ((heavy_amplicon < amplicons) &&
         (heavy_progress < heavy_amplicon_count))
    {
      unsigned int a = heavy_amplicon++;
      if (swarminfo[ampinfo[a].swarmid].mass >=
          static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++heavy_progress);
          pthread_mutex_unlock(&heavy_mutex);
          uint64_t m = 0;
          uint64_t v = 0;
          check_heavy_var(bloom_f, buffer1, a, &m, &v,
                          variant_list, variant_list2);
          pthread_mutex_lock(&heavy_mutex);
          heavy_variants += v;
        }
    }
  pthread_mutex_unlock(&heavy_mutex);
  xfree(buffer1);

  xfree(variant_list2);
  xfree(variant_list);
}

uint64_t mark_light_var(struct bloomflex_s * bloom,
                        unsigned int seed,
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
  uint64_t hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash, variant_list, & variant_count);

  for(auto i = 0U; i < variant_count; i++)
    bloomflex_set(bloom, variant_list[i].hash);

  return variant_count;
}

void mark_light_thread(int64_t t)
{
  (void) t;

  struct var_s * variant_list = static_cast<struct var_s *>
    (xmalloc(sizeof(struct var_s) * (7 * longestamplicon + 4)));

  pthread_mutex_lock(&light_mutex);
  while (light_progress < light_amplicon_count)
    {
      unsigned int a = light_amplicon--;
      if (swarminfo[ampinfo[a].swarmid].mass <
          static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++light_progress);
          pthread_mutex_unlock(&light_mutex);
          uint64_t v = mark_light_var(bloom_f, a, variant_list);
          pthread_mutex_lock(&light_mutex);
          light_variants += v;
        }
    }
  pthread_mutex_unlock(&light_mutex);

  xfree(variant_list);
}


/******************** FASTIDIOUS END ********************/


inline void find_variant_matches(unsigned int seed,
                                 var_s * var,
                                 unsigned int * hits_data,
                                 unsigned int * hits_count)
{
  if (! bloom_get(bloom_a, var->hash))
    return;

  /* compute hash and corresponding hash table index */

  uint64_t j = hash_getindex(var->hash);

  /* find matching buckets */

  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, var->hash))
        {
          unsigned int amp = hash_get_data(j);

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

void check_variants(unsigned int seed,
                    var_s * variant_list,
                    unsigned int * hits_data,
                    unsigned int * hits_count)
{
  unsigned int variant_count = 0;
  * hits_count = 0;

  char * sequence = db_getsequence(seed);
  unsigned int seqlen = db_getsequencelen(seed);
  uint64_t hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash, variant_list, & variant_count);

  for(auto i = 0U; i < variant_count; i++)
    find_variant_matches(seed, variant_list + i, hits_data, hits_count);
}

void network_thread(int64_t t)
{
  (void) t;

  unsigned int * hits_data = static_cast<unsigned int *>
    (xmalloc((7 * longestamplicon + 5) * sizeof(unsigned int)));

  struct var_s * variant_list = static_cast<struct var_s *>
    (xmalloc((7 * longestamplicon + 5) * sizeof(struct var_s)));

  pthread_mutex_lock(&network_mutex);
  while (network_amp < amplicons)
    {
      unsigned int amp = network_amp++;
      progress_update(amp);

      pthread_mutex_unlock(&network_mutex);

      unsigned int hits_count = 0;
      check_variants(amp, variant_list, hits_data, & hits_count);
      pthread_mutex_lock(&network_mutex);

      ampinfo[amp].link_start = network_count;
      ampinfo[amp].link_count = hits_count;

      if (network_count + hits_count > network_alloc)
        {
          while (network_count + hits_count > network_alloc)
            network_alloc += one_megabyte;
          network = static_cast<unsigned int*>
            (xrealloc(network, network_alloc * sizeof(unsigned int)));
        }

      for(auto i = 0U; i < hits_count; i++)
        network[network_count++] = hits_data[i];
    }
  pthread_mutex_unlock(&network_mutex);

  xfree(variant_list);
  xfree(hits_data);
}

void process_seed(unsigned int seed)
{
  /* update swarm stats */
  struct ampinfo_s * bp = ampinfo + seed;

  swarmsize++;
  if (bp->generation > swarm_maxgen)
    swarm_maxgen = bp->generation;
  uint64_t abundance = db_getabundance(seed);
  abundance_sum += abundance;
  if (abundance == 1)
    singletons++;
  swarm_sumlen += db_getsequencelen(seed);

  unsigned int s = ampinfo[seed].link_start;
  unsigned int c = ampinfo[seed].link_count;

  if (global_hits_count + c > global_hits_alloc)
    {
      while (global_hits_count + c > global_hits_alloc)
        global_hits_alloc += 4 * one_kilobyte;
      global_hits_data = static_cast<unsigned int *>
        (xrealloc(global_hits_data, global_hits_alloc * sizeof(unsigned int)));
    }

  for(auto i = 0U; i < c; i++)
    {
      unsigned int amp = network[s + i];

      if (ampinfo[amp].swarmid == no_swarm)
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
  /*
    Swarm checks that all amplicon sequences are unique (strictly
    dereplicated input data), so distinct amplicons with the same
    sequence are not expected at this stage.
  */

  const auto * x = static_cast<const unsigned int*>(a);
  const auto * y = static_cast<const unsigned int*>(b);
  int status {0};

  if (*x < *y) {
    status = -1;
  }
  else if (*x > *y) {
    status = +1;
  }
  else {
    status = 0;
  }
  return status;
}

int compare_mass(const void * a, const void * b)
{
  const swarminfo_s * x = swarminfo + *(static_cast<const unsigned int *>(a));
  const swarminfo_s * y = swarminfo + *(static_cast<const unsigned int *>(b));

  uint64_t m = x->mass;
  uint64_t n = y->mass;
  int status {0};

  if (m > n) {
    status = -1;
  }
  else if (m < n) {
    status = +1;
  }
  else {
    status = strcmp(db_getheader(x->seed), db_getheader(y->seed));
  }
  return status;
}

inline void add_amp_to_swarm(unsigned int amp)
{
  /* add to swarm */
  ampinfo[current_swarm_tail].next = amp;
  current_swarm_tail = amp;
}

void algo_d1_run()
{
  longestamplicon = db_getlongestsequence();
  amplicons = db_getsequencecount();

  ampinfo = static_cast<struct ampinfo_s *>
    (xmalloc(amplicons * sizeof(struct ampinfo_s)));

  global_hits_alloc = longestamplicon * 7 + 4 + 1;
  global_hits_data = static_cast<unsigned int *>
    (xmalloc(global_hits_alloc * sizeof(unsigned int)));

  /* compute hash for all amplicons and store them in a hash table */

  hash_alloc(amplicons);
  bloom_a = bloom_init(hash_get_tablesize());

  duplicates_found = 0;

  progress_init("Hashing sequences:", amplicons);
  for(auto i = 0U; i < amplicons; i++)
    {
      struct ampinfo_s * bp = ampinfo + i;
      bp->generation = 0;
      bp->swarmid = no_swarm;
      bp->next = no_swarm;
      bp->graft_cand = no_swarm;
      hash_insert(i);
      progress_update(i);
      if (duplicates_found)
        break;
    }

  if (duplicates_found)
    {
      fprintf(logfile,
              "\n\n"
              "Error: some fasta entries have identical sequences.\n"
              "Swarm expects dereplicated fasta files.\n"
              "Such files can be produced with swarm or vsearch:\n"
              " swarm -d 0 -w derep.fasta -o /dev/null input.fasta\n"
              "or\n"
              " vsearch --derep_fulllength input.fasta --sizein --sizeout --output derep.fasta\n");
      exit(1);
    }

  progress_done();

  unsigned char * dir = nullptr;
  uint64_t * hearray = nullptr;

  if (uclustfile)
    {
      dir = static_cast<unsigned char *>
        (xmalloc(longestamplicon*longestamplicon));
      hearray = static_cast<uint64_t *>
        (xmalloc(2 * longestamplicon * sizeof(uint64_t)));
    }

  /* for all amplicons, generate list of matching amplicons */

  network = static_cast<unsigned int*>
    (xmalloc(network_alloc * sizeof(unsigned int)));
  network_count = 0;

  pthread_mutex_init(&network_mutex, nullptr);
  network_amp = 0;
  progress_init("Building network: ", amplicons);
  ThreadRunner * network_tr = new ThreadRunner(static_cast<int>(opt_threads),
                                               network_thread);
  network_tr->run();
  delete network_tr;
  pthread_mutex_destroy(&network_mutex);

  progress_done();


  /* dump network to file */

  if (opt_network_file)
    {
      progress_init("Dumping network:  ", network_count);

      uint64_t n = 0;
      for(auto seed = 0U; seed < amplicons; seed++)
        {
          struct ampinfo_s * ap = ampinfo + seed;

          unsigned int link_start = ap->link_start;
          unsigned int link_count = ap->link_count;

          qsort(network + link_start,
                link_count,
                sizeof(unsigned int),
                compare_amp);

          for(auto i = 0U; i < link_count; i++)
            {
              unsigned int neighbour = network[link_start + i];
              fprint_id(network_file, seed);
              fprintf(network_file, "\t");
              fprint_id(network_file, neighbour);
              fprintf(network_file, "\n");
              n++;
            }
          progress_update(n);
        }
      progress_done();
    }


  /* for each non-swarmed amplicon look for subseeds ... */

  unsigned int swarmcount = 0;
  progress_init("Clustering:       ", amplicons);

  for(auto seed = 0U; seed < amplicons; seed++)
    {
      struct ampinfo_s * ap = ampinfo + seed;

      if (ap->swarmid == no_swarm)
        {
          /* start a new swarm with a new initial seed */

          ap->swarmid = swarmcount;
          ap->generation = 0;
          ap->parent = no_swarm;
          ap->next = no_swarm;

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
                sizeof(unsigned int), compare_amp);

          /* add subseeds on list to current swarm */
          for(auto i = 0U; i < global_hits_count; i++)
            add_amp_to_swarm(global_hits_data[i]);

          /* find later generation matches */
          unsigned int subseed = ap->next;
          while(subseed != no_swarm)
            {
              /* process all subseeds of this generation */
              global_hits_count = 0;

              while(subseed != no_swarm)
                {
                  process_seed(subseed);
                  subseed = ampinfo[subseed].next;
                }

              /* sort all of this generation */
              qsort(global_hits_data, global_hits_count,
                    sizeof(unsigned int), compare_amp);

              /* add them to the swarm */
              for(auto i = 0U; i < global_hits_count; i++)
                add_amp_to_swarm(global_hits_data[i]);

              /* start with most abundant amplicon of next generation */
              if (global_hits_count)
                subseed = global_hits_data[0];
              else
                subseed = no_swarm;
            }

          if (swarmcount >= swarminfo_alloc)
            {
              /* allocate memory for more swarms... */
              swarminfo_alloc += one_kilobyte;
              swarminfo = static_cast<struct swarminfo_s *>
                (xrealloc(swarminfo, swarminfo_alloc * sizeof(swarminfo_s)));
            }

          struct swarminfo_s * sp = swarminfo + swarmcount;

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

          swarmcount++;
        }
      progress_update(seed+1);
    }
  progress_done();

  xfree(global_hits_data);

  xfree(network);
  network = nullptr;

  swarmcount_adjusted = swarmcount;

  /* fastidious */

  if (opt_fastidious)
    {
      fprintf(logfile, "\n");
      fprintf(logfile, "Results before fastidious processing:\n");
      fprintf(logfile, "Number of swarms:  %u\n", swarmcount);
      fprintf(logfile, "Largest swarm:     %u\n", largest);
      fprintf(logfile, "\n");

      uint64_t small_otus = 0;
      uint64_t amplicons_in_small_otus = 0;
      uint64_t nucleotides_in_small_otus = 0;

      progress_init("Counting amplicons in heavy and light swarms",
                    swarmcount);

      for(auto i = 0ULL; i < swarmcount; i++)
        {
          struct swarminfo_s * sp = swarminfo + i;
          if (sp->mass < static_cast<uint64_t>(opt_boundary))
            {
              amplicons_in_small_otus += sp->size;
              nucleotides_in_small_otus += sp->sumlen;
              small_otus++;
            }
          progress_update(i+1);
        }
      progress_done();

      uint64_t amplicons_in_large_otus = amplicons - amplicons_in_small_otus;
      uint64_t large_otus = swarmcount - small_otus;

      fprintf(logfile, "Heavy swarms: %" PRIu64 ", with %" PRIu64 " amplicons\n",
              large_otus, amplicons_in_large_otus);
      fprintf(logfile, "Light swarms: %" PRIu64 ", with %" PRIu64 " amplicons\n",
              small_otus, amplicons_in_small_otus);
      fprintf(logfile, "Total length of amplicons in light swarms: %" PRIu64 "\n",
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

          uint64_t bits = static_cast<uint64_t>(opt_bloom_bits); /* 16 */

          // int64_t k = int(bits * 0.693);    /* 11 */
          unsigned int k = static_cast<unsigned int>(4 * bits / 10); /* 6 */
          if (k < 1)
            k = 1;

          uint64_t m = bits * 7 * nucleotides_in_small_otus;

          uint64_t memtotal = arch_get_memtotal();
          uint64_t memused = arch_get_memused();

          if (opt_ceiling)
            {
              uint64_t memrest
                = one_megabyte * static_cast<uint64_t>(opt_ceiling) - memused;
              uint64_t new_bits = 8 * memrest / (7 * nucleotides_in_small_otus);
              if (new_bits < bits)
                {
                  if (new_bits < 2)
                    fatal("Insufficient memory remaining for Bloom filter");
                  fprintf(logfile, "Reducing memory used for Bloom filter due to --ceiling option.\n");
                  bits = new_bits;
                  // k = int(bits * 0.693);
                  k = static_cast<unsigned int>(4 * bits / 10);
                  if (k < 1)
                    k = 1;

                  m = bits * 7 * nucleotides_in_small_otus;
                }
            }

          if (m < 64)
            m = 64;

          if (memused + m/8 > memtotal)
            {
              fprintf(logfile, "WARNING: Memory usage will probably exceed total amount of memory available.\n");
              fprintf(logfile, "Try to reduce memory footprint using the --bloom-bits or --ceiling options.\n");
            }

          fprintf(logfile,
                  "Bloom filter: bits=%" PRIu64 ", m=%" PRIu64 ", k=%u, size=%.1fMB\n",
                  bits, m, k, static_cast<double>(m) / (8 * one_megabyte));


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

          pthread_mutex_init(&light_mutex, nullptr);
          light_progress = 0;
          light_amplicon_count = amplicons_in_small_otus;
          light_amplicon = amplicons - 1;
          ThreadRunner * tr = new ThreadRunner(static_cast<int>(opt_threads),
                                               mark_light_thread);
          tr->run();
          delete tr;
          pthread_mutex_destroy(&light_mutex);

          progress_done();

          fprintf(logfile,
                  "Generated %" PRIu64 " variants from light swarms\n",
                  light_variants);

          progress_init("Checking heavy swarm amplicons against Bloom filter",
                        amplicons_in_large_otus);

          /* process amplicons in order from most to least abundant */
          /* but stop when all amplicons in large otus are processed */

          pthread_mutex_init(&graft_mutex, nullptr);

          heavy_variants = 0;

          pthread_mutex_init(&heavy_mutex, nullptr);
          heavy_progress = 0;
          heavy_amplicon_count = amplicons_in_large_otus;
          heavy_amplicon = 0;
          ThreadRunner * heavy_tr
            = new ThreadRunner(static_cast<int> (opt_threads),
                               check_heavy_thread);
          heavy_tr->run();
          delete heavy_tr;
          pthread_mutex_destroy(&heavy_mutex);

          progress_done();

          bloomflex_exit(bloom_f);

          pthread_mutex_destroy(&graft_mutex);

          fprintf(logfile, "Heavy variants: %" PRIu64 "\n", heavy_variants);
          fprintf(logfile, "Got %" PRId64 " graft candidates\n", graft_candidates);
          unsigned int grafts = attach_candidates(amplicons);
          fprintf(logfile, "Made %u grafts\n", grafts);
          fprintf(logfile, "\n");
        }
    }


  /* dump swarms */

  progress_init("Writing swarms:   ", swarmcount);

  if (opt_mothur)
    fprintf(outfile, "swarm_%" PRId64 "\t%" PRIu64, opt_differences, swarmcount_adjusted);

  for(auto i = 0U; i < swarmcount; i++)
    {
      if (!swarminfo[i].attached)
        {
          unsigned int seed = swarminfo[i].seed;
          for(auto a = seed; a != no_swarm; a = ampinfo[a].next)
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
                    fputc(sepchar, outfile);
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

      unsigned int * sorter = static_cast<unsigned int *>
        (xmalloc(swarmcount * sizeof(unsigned int)));
      for(auto i = 0U; i < swarmcount; i++)
        sorter[i] = i;
      qsort(sorter, swarmcount, sizeof(unsigned int), compare_mass);

      for(auto j = 0U; j < swarmcount; j++)
        {
          unsigned int i = sorter[j];
          if (!swarminfo[i].attached)
            {
              unsigned int seed = swarminfo[i].seed;
              fprintf(fp_seeds, ">");
              fprint_id_with_new_abundance(fp_seeds, seed, swarminfo[i].mass);
              fprintf(fp_seeds, "\n");
              db_fprintseq(fp_seeds, seed, 0);
            }
          progress_update(i+1);
        }

      xfree(sorter);

      progress_done();
    }


  /* output internal structure */

  if (opt_internal_structure)
    {
      unsigned int cluster_no = 0;

      progress_init("Writing structure:", swarmcount);

      for(auto swarmid = 0U; swarmid < swarmcount; swarmid++)
        {
          if (!swarminfo[swarmid].attached)
            {
              unsigned int seed = swarminfo[swarmid].seed;

              struct ampinfo_s * bp = ampinfo + seed;

              for(auto a = bp->next; a != no_swarm; a = ampinfo[a].next)
                {
                  uint64_t graft_parent = ampinfo[a].graft_cand;
                  if (graft_parent != no_swarm)
                    {
                      fprint_id_noabundance(internal_structure_file,
                                            graft_parent);
                      fprintf(internal_structure_file, "\t");
                      fprint_id_noabundance(internal_structure_file, a);
                      fprintf(internal_structure_file,
                              "\t%d\t%u\t%u\n",
                              2,
                              cluster_no + 1,
                              ampinfo[graft_parent].generation + 1);
                    }

                  uint64_t parent = ampinfo[a].parent;
                  if (parent != no_swarm)
                    {
                      fprint_id_noabundance(internal_structure_file, parent);
                      fprintf(internal_structure_file, "\t");
                      fprint_id_noabundance(internal_structure_file, a);
                      fprintf(internal_structure_file,
                              "\t%u\t%u\t%u\n",
                              1U,
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

      for(auto swarmid = 0U; swarmid < swarmcount; swarmid++)
        {
          if (!swarminfo[swarmid].attached)
            {
              unsigned int seed = swarminfo[swarmid].seed;

              struct ampinfo_s * bp = ampinfo + seed;

              fprintf(uclustfile, "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                      cluster_no,
                      swarminfo[swarmid].size);
              fprint_id(uclustfile, seed);
              fprintf(uclustfile, "\t*\n");

              fprintf(uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                      cluster_no,
                      db_getsequencelen(seed));
              fprint_id(uclustfile, seed);
              fprintf(uclustfile, "\t*\n");

              for(auto a = bp->next; a != no_swarm; a = ampinfo[a].next)
                {
                  char * dseq = db_getsequence(a);
                  int64_t dlen = db_getsequencelen(a);
                  char * qseq = db_getsequence(seed);
                  int64_t qlen = db_getsequencelen(seed);

                  int64_t nwscore = 0;
                  int64_t nwdiff = 0;
                  char * nwalignment = nullptr;
                  int64_t nwalignmentlength = 0;

                  nw(dseq, dlen, qseq, qlen,
                     score_matrix_63, penalty_gapopen, penalty_gapextend,
                     & nwscore, & nwdiff, & nwalignmentlength, & nwalignment,
                     dir, reinterpret_cast<int64_t *>(hearray), 0, 0);

                  double percentid
                    = 100.0 * static_cast<double>(nwalignmentlength - nwdiff)
                    / static_cast<double>(nwalignmentlength);

                  fprintf(uclustfile,
                          "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                          cluster_no,
                          db_getsequencelen(a),
                          percentid,
                          nwdiff > 0 ? nwalignment : "=");

                  fprint_id(uclustfile, a);
                  fprintf(uclustfile, "\t");
                  fprint_id(uclustfile, seed);
                  fprintf(uclustfile, "\n");

                  if (nwalignment)
                    xfree(nwalignment);
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
      for(auto i = 0ULL; i < swarmcount; i++)
        {
          swarminfo_s * sp = swarminfo + i;
          if (!sp->attached)
            {
              fprintf(statsfile, "%u\t%" PRIu64 "\t", sp->size, sp->mass);
              fprint_id_noabundance(statsfile, sp->seed);
              fprintf(statsfile, "\t%" PRIu64 "\t%u\t%u\t%u\n",
                      db_getabundance(sp->seed),
                      sp->singletons, sp->maxgen, sp->maxgen);
            }
          progress_update(i);
        }
      progress_done();
    }

  fprintf(logfile, "\n");
  fprintf(logfile, "Number of swarms:  %" PRIu64 "\n", swarmcount_adjusted);
  fprintf(logfile, "Largest swarm:     %u\n", largest);
  fprintf(logfile, "Max generations:   %u\n", maxgen);

  bloom_exit(bloom_a);
  hash_free();

  if (swarminfo)
    xfree(swarminfo);

  xfree(ampinfo);

  if (uclustfile)
    {
      xfree(dir);
      xfree(hearray);
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
