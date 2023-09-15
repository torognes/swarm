/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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
#include "arch.h"
#include "bloomflex.h"
#include "bloompat.h"
#include "db.h"
#include "hashtable.h"
#include "nw.h"
#include "variants.h"
#include "util.h"
#include "utils/nt_codec.h"
#include "utils/score_matrix.h"
#include "utils/threads.h"
#include "zobrist.h"
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <climits>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc(), size_t
#include <cstdlib>  // qsort()
#include <cstring>  // std::memcmp
#include <pthread.h>
#include <vector>


#ifndef PRIu64
#ifdef _WIN32
#define PRIu64 "I64u"
#else
constexpr char PRIu64[] = "lu";
#endif
#endif

#ifndef PRId64
#ifdef _WIN32
#define PRId64 "I64d"
#else
constexpr char PRId64[] = "ld";
#endif
#endif


constexpr unsigned int one_kilobyte {1U << 10U};
constexpr unsigned int one_megabyte {one_kilobyte * one_kilobyte};

static uint64_t duplicates_found {0};  // several function calls

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

/* Information about each swarm (cluster) */

static struct swarminfo_s
{
  uint64_t mass {0}; /* the sum of abundances of amplicons in this swarm */
  uint64_t sumlen {0}; /* sum of length of amplicons in swarm */
  unsigned int seed {0}; /* amplicon id of the initial seed of this swarm */
  unsigned int last {0}; /* amplicon id of the last seed in this swarm */
  unsigned int size {0}; /* total number of amplicons in this swarm */
  unsigned int singletons {0}; /* number of amplicons with abundance 1 */
  unsigned int maxgen {0}; /* the generation of the amplicon farthest from seed */
  bool attached {false}; /* this is a small swarm attached to a large (fastidious) */
  char dummy_1 = '\0'; /* alignment padding only */
  char dummy_2 = '\0'; /* alignment padding only */
  char dummy_3 = '\0'; /* alignment padding only */
} * swarminfo = nullptr;

static struct graft_cand
{
  unsigned int parent;
  unsigned int child;
} * graft_array = nullptr;

/* Information about potential grafts */
static int64_t graft_candidates {0};
static pthread_mutex_t graft_mutex;

constexpr unsigned int no_swarm {UINT_MAX};

static unsigned int current_swarm_tail {0};

static uint64_t swarminfo_alloc {0};

/* overall statistics */
static unsigned int maxgen {0};
static unsigned int largest {0};
static uint64_t swarmcount_adjusted {0};

/* per swarm statistics */
static unsigned int singletons {0};
static uint64_t abundance_sum {0}; /* = mass */
static unsigned int swarmsize {0};
static unsigned int swarm_maxgen {0};
static uint64_t swarm_sumlen {0};

static unsigned int * global_hits_data {nullptr};
static unsigned int global_hits_alloc {0};
static unsigned int global_hits_count {0};

static unsigned int longestamplicon {0};

static unsigned int amplicons {0};

static pthread_mutex_t heavy_mutex;
static uint64_t heavy_variants {0};
static uint64_t heavy_progress {0};
static uint64_t heavy_amplicon_count {0};
static unsigned int heavy_amplicon {0};

static pthread_mutex_t light_mutex;
static uint64_t light_variants {0};
static uint64_t light_progress {0};
static uint64_t light_amplicon_count {0};
static unsigned int light_amplicon {0};

static uint64_t network_alloc {one_megabyte};
static unsigned int * network {nullptr};
static unsigned int network_count {0};
static pthread_mutex_t network_mutex;
static unsigned int network_amp {0};

static struct bloom_s * bloom_a {nullptr}; // Bloom filter for amplicons

static struct bloomflex_s * bloom_f {nullptr}; // Huge Bloom filter for fastidious

void attach(unsigned int seed, unsigned int amp);
void add_graft_candidate(unsigned int seed, unsigned int amp);
auto compare_grafts(const void * a, const void * b) -> int;
auto attach_candidates(unsigned int amplicon_count) -> unsigned int;
auto hash_check_attach(char * seed_sequence,
                       unsigned int seed_seqlen,
                       struct var_s * var,
                       unsigned int seed) -> bool;
void check_heavy_var(struct bloomflex_s * bloom,
                     char * varseq,
                     unsigned int seed,
                     uint64_t * m,
                     uint64_t * v,
                     struct var_s * variant_list,
                     struct var_s * variant_list2);
void check_heavy_thread(int64_t t);
auto mark_light_var(struct bloomflex_s * bloom,
                        unsigned int seed,
                        struct var_s * variant_list) -> uint64_t;
void mark_light_thread(int64_t t);
void check_variants(unsigned int seed,
                    var_s * variant_list,
                    unsigned int * hits_data,
                    unsigned int * hits_count);
void network_thread(int64_t t);
void process_seed(unsigned int seed);
auto compare_amp(const void * a, const void * b) -> int;
auto compare_mass(const void * a, const void * b) -> int;


inline auto check_amp_identical(unsigned int amp1,
                                unsigned int amp2) -> bool
{
  /* amplicon are identical if they have the same length, and the
     exact same sequence */
  const unsigned int amp1_seqlen = db_getsequencelen(amp1);
  const unsigned int amp2_seqlen = db_getsequencelen(amp2);

  return ((amp1_seqlen == amp2_seqlen) and
          (std::memcmp(db_getsequence(amp1),
                  db_getsequence(amp2),
                  nt_bytelength(amp1_seqlen)) == 0));
}


inline void hash_insert(unsigned int amp)
{
  /* find the first empty bucket */
  const uint64_t hash = db_gethash(amp);
  uint64_t j = hash_getindex(hash);
  bool duplicate {false};
  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, hash) and
          check_amp_identical(amp, hash_get_data(j))) {
        duplicate = true;
      }
      j = hash_getnextindex(j);
    }

  if (duplicate) {
    ++duplicates_found;
  }

  hash_set_occupied(j);
  hash_set_value(j, hash);
  hash_set_data(j, amp);

  bloom_set(bloom_a, hash);
}


/******************** FASTIDIOUS START ********************/


auto attach(unsigned int seed, unsigned int amp) -> void
{
  /* graft light swarm (amp) on heavy swarm (seed) */

  swarminfo_s * heavy_swarm = swarminfo + ampinfo[seed].swarmid;
  swarminfo_s * light_swarm = swarminfo + ampinfo[amp].swarmid;

  // attach the seed of the light swarm to the tail of the heavy swarm (refactoring: unclear)
  ampinfo[heavy_swarm->last].next = light_swarm->seed;
  heavy_swarm->last = light_swarm->last;

  // Update swarm info
  heavy_swarm->size += light_swarm->size;
  heavy_swarm->singletons += light_swarm->singletons;
  heavy_swarm->mass += light_swarm->mass;
  heavy_swarm->sumlen += light_swarm->sumlen;
  /* maxgen is untouched */

  /* flag attachment to avoid doing it again */
  light_swarm->attached = true;

  // Update overall stats
  if (heavy_swarm->size > largest) {
    largest = heavy_swarm->size;
  }

  --swarmcount_adjusted;
}


auto add_graft_candidate(unsigned int seed, unsigned int amp) -> void
{
  pthread_mutex_lock(&graft_mutex);
  ++graft_candidates;
  // if there is no heavy candidate to graft amp, or if seed is
  // earlier in the sorting order, then we change the attachment to
  // seed
  if ((ampinfo[amp].graft_cand == no_swarm) or (ampinfo[amp].graft_cand > seed)) {
    ampinfo[amp].graft_cand = seed;
  }
  pthread_mutex_unlock(&graft_mutex);
}


auto compare_grafts(const void * a, const void * b) -> int
{
  const auto * x = static_cast<const struct graft_cand *>(a);
  const auto * y = static_cast<const struct graft_cand *>(b);
  int status {0};

  if (x->parent < y->parent) {
    status = -1;
  }
  else if (x->parent > y->parent) {
    status = +1;
  }
  else { // x->parent == y->parent
    // assert(x->child != y->child); always true, children have different index values, by construction
    // assert(x->child < y->child);  that seems to be always true in my tests, even on large datasets
    if (x->child < y->child) {
      status = -1;
    }
    else if (x->child > y->child) {  // unreachable?
      status = +1;
    }
    else {
      status = 0;  // unreachable, replace with above assertion
    }
  }

  return status;
}


auto attach_candidates(unsigned int amplicon_count) -> unsigned int
{
  /* count pairs */
  unsigned int pair_count = 0;
  for(auto i = 0U; i < amplicon_count; i++) {
    if (ampinfo[i].graft_cand != no_swarm) {
      ++pair_count;
    }
  }

  unsigned int grafts = 0;
  progress_init("Grafting light swarms on heavy swarms", pair_count);

  /* allocate memory */
  graft_array = new struct graft_cand[pair_count];  // refactor to std::vector

  /* fill in */
  unsigned int j = 0;
  for(auto i = 0U; i < amplicon_count; i++) {
    if (ampinfo[i].graft_cand != no_swarm)
      {
        graft_array[j].parent = ampinfo[i].graft_cand;
        graft_array[j].child = i;  // so two children cannot have the same uint value
        ++j;
      }
  }

  /* sort */
  std::qsort(graft_array, pair_count, sizeof(struct graft_cand), compare_grafts);

  /* attach in order */
  for(auto i = 0U; i < pair_count; i++)
    {
      const unsigned int parent = graft_array[i].parent;
      const unsigned int child  = graft_array[i].child;

      if (swarminfo[ampinfo[child].swarmid].attached)
        {
          /* this light swarm is already attached */
          ampinfo[child].graft_cand = no_swarm;
        }
      else
        {
          /* attach child to parent */
          attach(parent, child);
          ++grafts;
        }
      progress_update(i+1);
    }
  progress_done();
  delete [] graft_array;
  graft_array = nullptr;
  return grafts;
}


auto hash_check_attach(char * seed_sequence,
                       unsigned int seed_seqlen,
                       struct var_s * var,
                       unsigned int seed) -> bool
{
  /* seed is the original large swarm seed */

  /* compute hash and corresponding hash table index */
  const uint64_t hash = var->hash;
  uint64_t j = hash_getindex(hash);

  /* find matching buckets */

  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, hash))
        {
          /* check that mass is below threshold */
          const unsigned int amp = hash_get_data(j);

          /* make absolutely sure sequences are identical */
          char * amp_sequence = db_getsequence(amp);
          const unsigned int amp_seqlen = db_getsequencelen(amp);
          if (check_variant(seed_sequence, seed_seqlen, var, amp_sequence, amp_seqlen))
            {
              add_graft_candidate(seed, amp);
              return true;
            }
        }
      j = hash_getnextindex(j);
    }
  return false;
}


inline auto check_heavy_var_2(char * seq,
                              unsigned int seqlen,
                              unsigned int seed,
                              struct var_s * variant_list) -> uint64_t
{
  /* Check second generation microvariants of the heavy swarm amplicons
     and see if any of them are identical to a light swarm amplicon. */

  uint64_t matches = 0;
  unsigned int variant_count = 0;

  const uint64_t hash = zobrist_hash(reinterpret_cast<unsigned char *>(seq), seqlen);
  generate_variants(seq, seqlen, hash, variant_list, &variant_count);

  for(auto i = 0U; i < variant_count; i++) {
    if (bloom_get(bloom_a, variant_list[i].hash) and
        hash_check_attach(seq, seqlen, variant_list + i, seed)) {
      ++matches;
    }
  }

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
  const unsigned int seqlen = db_getsequencelen(seed);
  const uint64_t hash = db_gethash(seed);
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
  static constexpr unsigned int i {7};  // max number of microvariants = 7 * len + 4
  static constexpr unsigned int j {4};  //                               i * len + j
  static constexpr unsigned int nt_per_uint64 {32};  // 32 nucleotides can fit in a uint64
  (void) t;

  auto * variant_list = new struct var_s[i * longestamplicon + j];
  auto * variant_list2 = new struct var_s[i * (longestamplicon + 1) + j];

  const std::size_t size =
    sizeof(uint64_t) * ((db_getlongestsequence() + 2 + nt_per_uint64 - 1) / nt_per_uint64);
  std::vector<char> buffer1(size);
  pthread_mutex_lock(&heavy_mutex);
  while ((heavy_amplicon < amplicons) and
         (heavy_progress < heavy_amplicon_count))
    {
      const unsigned int heavy_amplicon_id = heavy_amplicon++;
      if (swarminfo[ampinfo[heavy_amplicon_id].swarmid].mass >=
          static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++heavy_progress);
          pthread_mutex_unlock(&heavy_mutex);
          uint64_t m {0};
          uint64_t v {0};
          check_heavy_var(bloom_f, buffer1.data(), heavy_amplicon_id, &m, &v,
                          variant_list, variant_list2);
          pthread_mutex_lock(&heavy_mutex);
          heavy_variants += v;
        }
    }
  pthread_mutex_unlock(&heavy_mutex);
  delete [] variant_list2;
  variant_list2 = nullptr;
  delete [] variant_list;
  variant_list = nullptr;
}


auto mark_light_var(struct bloomflex_s * bloom,
                        unsigned int seed,
                        struct var_s * variant_list) -> uint64_t
{
  /*
    add all microvariants of seed to Bloom filter

    bloom is a BloomFilter in which to enter the variants
    seed is the original seed
  */

  hash_insert(seed);

  unsigned int variant_count = 0;

  char * sequence = db_getsequence(seed);
  const unsigned int seqlen = db_getsequencelen(seed);
  const uint64_t hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash, variant_list, & variant_count);

  for(auto i = 0U; i < variant_count; i++) {
    bloomflex_set(bloom, variant_list[i].hash);
  }

  return variant_count;
}


void mark_light_thread(int64_t t)
{
  static constexpr unsigned int i {7};  // max number of microvariants = 7 * len + 4
  static constexpr unsigned int j {4};  //                               i * len + j

  (void) t;

  auto * variant_list = new struct var_s[i * longestamplicon + j];

  pthread_mutex_lock(&light_mutex);
  while (light_progress < light_amplicon_count)
    {
      const unsigned int light_amplicon_id = light_amplicon--;
      if (swarminfo[ampinfo[light_amplicon_id].swarmid].mass <
          static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++light_progress);
          pthread_mutex_unlock(&light_mutex);
          const uint64_t v = mark_light_var(bloom_f, light_amplicon_id, variant_list);
          pthread_mutex_lock(&light_mutex);
          light_variants += v;
        }
    }
  pthread_mutex_unlock(&light_mutex);

  delete [] variant_list;
  variant_list = nullptr;
}


/******************** FASTIDIOUS END ********************/


inline void find_variant_matches(unsigned int seed,
                                 var_s * var,
                                 unsigned int * hits_data,
                                 unsigned int * hits_count)
{
  if (not bloom_get(bloom_a, var->hash)) {
    return;
  }

  /* compute hash and corresponding hash table index */

  uint64_t j = hash_getindex(var->hash);

  /* find matching buckets */

  while (hash_is_occupied(j))
    {
      if (hash_compare_value(j, var->hash))
        {
          const unsigned int amp = hash_get_data(j);

          /* avoid self */
          if (seed != amp) {
            if ((opt_no_otu_breaking) or
                (db_getabundance(seed) >= db_getabundance(amp)))
              {
                char * seed_sequence = db_getsequence(seed);
                const unsigned int seed_seqlen = db_getsequencelen(seed);

                char * amp_sequence = db_getsequence(amp);
                const unsigned int amp_seqlen = db_getsequencelen(amp);

                if (check_variant(seed_sequence, seed_seqlen,
                                  var,
                                  amp_sequence, amp_seqlen))
                  {
                    hits_data[(*hits_count)++] = amp;
                    break;
                  }
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
  const unsigned int seqlen = db_getsequencelen(seed);
  const uint64_t hash = db_gethash(seed);
  generate_variants(sequence, seqlen, hash, variant_list, & variant_count);

  for(auto i = 0U; i < variant_count; i++) {
    find_variant_matches(seed, variant_list + i, hits_data, hits_count);
  }
}


void network_thread(int64_t t)
{
  static constexpr unsigned int i {7};  // max number of microvariants = 7 * len + 4
  static constexpr unsigned int j {4};  //                               i * len + j

  (void) t;

  auto * hits_data = new unsigned int[i * longestamplicon + j + 1];
  auto * variant_list = new struct var_s[i * longestamplicon + j + 1];

  pthread_mutex_lock(&network_mutex);
  while (network_amp < amplicons)
    {
      const unsigned int amp = network_amp++;
      progress_update(amp);

      pthread_mutex_unlock(&network_mutex);

      unsigned int hits_count = 0;
      check_variants(amp, variant_list, hits_data, & hits_count);
      pthread_mutex_lock(&network_mutex);

      ampinfo[amp].link_start = network_count;
      ampinfo[amp].link_count = hits_count;

      if (network_count + hits_count > network_alloc)
        {
          while (network_count + hits_count > network_alloc) {
            network_alloc += one_megabyte;
          }
          network = static_cast<unsigned int*>
            (xrealloc(network, network_alloc * sizeof(unsigned int)));
        }

      for(auto k = 0U; k < hits_count; k++) {
        network[network_count++] = hits_data[k];
      }
    }
  pthread_mutex_unlock(&network_mutex);

  delete [] variant_list;
  variant_list = nullptr;
  delete [] hits_data;
  hits_data = nullptr;
}


void process_seed(unsigned int seed)
{
  /* update swarm stats */
  struct ampinfo_s * bp = ampinfo + seed;

  ++swarmsize;
  if (bp->generation > swarm_maxgen) {
    swarm_maxgen = bp->generation;
  }
  const uint64_t abundance = db_getabundance(seed);
  abundance_sum += abundance;
  if (abundance == 1) {
    ++singletons;
  }
  swarm_sumlen += db_getsequencelen(seed);

  const unsigned int link_start = ampinfo[seed].link_start;
  const unsigned int link_count = ampinfo[seed].link_count;

  if (global_hits_count + link_count > global_hits_alloc)
    {
      while (global_hits_count + link_count > global_hits_alloc) {
        global_hits_alloc += 4 * one_kilobyte;
      }
      global_hits_data = static_cast<unsigned int *>
        (xrealloc(global_hits_data, global_hits_alloc * sizeof(unsigned int)));
    }

  for(auto offset = 0U; offset < link_count; offset++)
    {
      const unsigned int amp = network[link_start + offset];

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


auto compare_amp(const void * a, const void * b) -> int
{
  /*
    earlier steps in swarm check that all amplicon sequences are
    unique (strictly dereplicated input data), so distinct amplicons
    with the same sequence are not expected at this stage.
  */
  const auto * x = static_cast<const unsigned int*>(a);
  const auto * y = static_cast<const unsigned int*>(b);
  int status {1};  // default is *x > *y

  assert(*x != *y);  // '*x == *y' is not expected at that stage

  if (*x < *y) {
    status = -1;
  }

  return status;
}


auto compare_mass(const void * a, const void * b) -> int
{
  const swarminfo_s * x = swarminfo + *(static_cast<const unsigned int *>(a));
  const swarminfo_s * y = swarminfo + *(static_cast<const unsigned int *>(b));

  const uint64_t m = x->mass;
  const uint64_t n = y->mass;
  int status {0};

  if (m > n) {
    status = -1;
  }
  else if (m < n) {
    status = +1;
  }
  else {
    status = std::strcmp(db_getheader(x->seed), db_getheader(y->seed));
  }
  return status;
}


inline void add_amp_to_swarm(unsigned int amp)
{
  /* add to swarm */
  ampinfo[current_swarm_tail].next = amp;
  current_swarm_tail = amp;
}


auto write_network_file(const unsigned int network_count,
                        struct Parameters const & parameters) -> void {
  progress_init("Dumping network:  ", network_count);

  uint64_t n_processed = 0;  // refactoring: reduce scope (add to for loop init)
  for(auto seed = 0U; seed < amplicons; seed++)
    {
      struct ampinfo_s * ap = ampinfo + seed;

      const unsigned int link_start = ap->link_start;
      const unsigned int link_count = ap->link_count;

      std::qsort(network + link_start,
            link_count,
            sizeof(unsigned int),
            compare_amp);

      for(auto link = 0U; link < link_count; link++)
        {
          const unsigned int neighbour = network[link_start + link];
          fprint_id(network_file, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(network_file, "\t");
          fprint_id(network_file, neighbour, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(network_file, "\n");
          ++n_processed;
        }
      progress_update(n_processed);
    }
  progress_done();
}


auto write_swarms_default_format(const unsigned int swarmcount,
                                 struct Parameters const & parameters) -> void {
  progress_init("Writing swarms:   ", swarmcount);

  for(auto i = 0U; i < swarmcount; i++) {
    if (not swarminfo[i].attached) {
      const unsigned int seed = swarminfo[i].seed;
      for(auto a = seed; a != no_swarm; a = ampinfo[a].next) {
        if (a != seed) {
          std::fputc(sepchar, outfile);
        }
        fprint_id(outfile, a,
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      }
      std::fputc('\n', outfile);
    }
    progress_update(i + 1);
  }

  progress_done();
}


auto write_swarms_mothur_format(const unsigned int swarmcount,
                                struct Parameters const & parameters) -> void {
  progress_init("Writing swarms:   ", swarmcount);

  std::fprintf(outfile, "swarm_%" PRId64 "\t%" PRIu64,
          parameters.opt_differences, swarmcount_adjusted);

  for(auto i = 0U; i < swarmcount; i++) {
    if (not swarminfo[i].attached) {
      const unsigned int seed = swarminfo[i].seed;
      for(auto a = seed; a != no_swarm; a = ampinfo[a].next) {
        if (a == seed) {
          std::fputc('\t', outfile);
        }
        else {
          std::fputc(',', outfile);
        }
        fprint_id(outfile, a,
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      }
    }
    progress_update(i + 1);
  }

  std::fputc('\n', outfile);

  progress_done();
}


auto write_swarms_uclust_format(const unsigned int swarmcount,
                                struct Parameters const & parameters,
                                unsigned char * dir,
                                uint64_t * hearray) -> void {
  static constexpr unsigned int one_hundred {100};
  unsigned int cluster_no = 0;
  const auto score_matrix_63 = create_score_matrix<int64_t>(parameters.penalty_mismatch);
  dir = new unsigned char[longestamplicon * longestamplicon];
  hearray = new uint64_t[2 * longestamplicon];

  progress_init("Writing UCLUST:   ", swarmcount);

  for(auto swarmid = 0U; swarmid < swarmcount; swarmid++)
    {
      if (not swarminfo[swarmid].attached)
        {
          const unsigned int seed = swarminfo[swarmid].seed;

          struct ampinfo_s * bp = ampinfo + seed;

          std::fprintf(uclustfile, "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                  cluster_no,
                  swarminfo[swarmid].size);
          fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(uclustfile, "\t*\n");

          std::fprintf(uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                  cluster_no,
                  db_getsequencelen(seed));
          fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(uclustfile, "\t*\n");

          for(auto a = bp->next; a != no_swarm; a = ampinfo[a].next)
            {
              char * dseq = db_getsequence(a);
              const int64_t dlen = db_getsequencelen(a);
              char * qseq = db_getsequence(seed);
              const int64_t qlen = db_getsequencelen(seed);

              int64_t nwscore = 0;
              int64_t nwdiff = 0;
              char * nwalignment = nullptr;
              int64_t nwalignmentlength = 0;

              nw(dseq, dlen, qseq, qlen,
                 score_matrix_63, penalty_gapopen, penalty_gapextend,
                 & nwscore, & nwdiff, & nwalignmentlength, & nwalignment,
                 dir, reinterpret_cast<int64_t *>(hearray), 0, 0);

              const double percentid
                = one_hundred * static_cast<double>(nwalignmentlength - nwdiff)
                / static_cast<double>(nwalignmentlength);

              std::fprintf(uclustfile,
                      "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                      cluster_no,
                      db_getsequencelen(a),
                      percentid,
                      nwdiff > 0 ? nwalignment : "=");

              fprint_id(uclustfile, a, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
              std::fprintf(uclustfile, "\t");
              fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
              std::fprintf(uclustfile, "\n");

              if (nwalignment != nullptr) {
                xfree(nwalignment);
              }
              nwalignment = nullptr;
            }

          ++cluster_no;
        }
      progress_update(swarmid);
    }
  progress_done();
  delete [] dir;
  dir = nullptr;
  delete [] hearray;
  hearray = nullptr;
}


auto write_representative_sequences(const unsigned int swarmcount,
                                    struct Parameters const & parameters) -> void {
  progress_init("Writing seeds:    ", swarmcount);

  auto * sorter = new unsigned int[swarmcount];
  for(auto i = 0U; i < swarmcount; i++) {
    sorter[i] = i;
  }
  std::qsort(sorter, swarmcount, sizeof(unsigned int), compare_mass);

  for(auto j = 0U; j < swarmcount; j++)
    {
      const unsigned int i = sorter[j];
      if (not swarminfo[i].attached)
        {
          const unsigned int seed = swarminfo[i].seed;
          std::fprintf(fp_seeds, ">");
          fprint_id_with_new_abundance(fp_seeds, seed, swarminfo[i].mass, parameters.opt_usearch_abundance);
          std::fprintf(fp_seeds, "\n");
          db_fprintseq(fp_seeds, seed);
        }
      progress_update(i + 1);
    }

  delete [] sorter;
  sorter = nullptr;

  progress_done();
}


auto write_structure_file(const unsigned int swarmcount,
                          struct Parameters const & parameters) -> void {
  unsigned int cluster_no {0};

  progress_init("Writing structure:", swarmcount);

  for(auto swarmid = 0U; swarmid < swarmcount; swarmid++)
    {
      if (not swarminfo[swarmid].attached)
        {
          const unsigned int seed = swarminfo[swarmid].seed;

          struct ampinfo_s * bp = ampinfo + seed;

          for(auto a = bp->next; a != no_swarm; a = ampinfo[a].next)
            {
              const uint64_t graft_parent = ampinfo[a].graft_cand;
              if (graft_parent != no_swarm)
                {
                  fprint_id_noabundance(internal_structure_file,
                                        graft_parent, parameters.opt_usearch_abundance);
                  std::fprintf(internal_structure_file, "\t");
                  fprint_id_noabundance(internal_structure_file, a, parameters.opt_usearch_abundance);
                  std::fprintf(internal_structure_file,
                          "\t%d\t%u\t%u\n",
                          2,
                          cluster_no + 1,
                          ampinfo[graft_parent].generation + 1);
                }

              const uint64_t parent = ampinfo[a].parent;
              if (parent != no_swarm)
                {
                  fprint_id_noabundance(internal_structure_file, parent, parameters.opt_usearch_abundance);
                  std::fprintf(internal_structure_file, "\t");
                  fprint_id_noabundance(internal_structure_file, a, parameters.opt_usearch_abundance);
                  std::fprintf(internal_structure_file,
                          "\t%u\t%u\t%u\n",
                          1U,
                          cluster_no + 1,
                          ampinfo[a].generation);
                }
            }

          ++cluster_no;
        }
      progress_update(swarmid);
    }
  progress_done();
}


auto write_stats_file(const unsigned int swarmcount,
                      struct Parameters const & parameters) -> void {
  progress_init("Writing stats:    ", swarmcount);
  for(auto i = 0ULL; i < swarmcount; i++)
    {
      swarminfo_s * sp = swarminfo + i;
      if (not sp->attached)
        {
          std::fprintf(statsfile, "%u\t%" PRIu64 "\t", sp->size, sp->mass);
          fprint_id_noabundance(statsfile, sp->seed, parameters.opt_usearch_abundance);
          std::fprintf(statsfile, "\t%" PRIu64 "\t%u\t%u\t%u\n",
                  db_getabundance(sp->seed),
                  sp->singletons, sp->maxgen, sp->maxgen);
        }
      progress_update(i);
    }
  progress_done();
}


void algo_d1_run(struct Parameters const & parameters)
{
  longestamplicon = db_getlongestsequence();
  amplicons = db_getsequencecount();

  ampinfo = new struct ampinfo_s[amplicons];

  // max number of microvariants = 7 * len + 4
  static constexpr unsigned int m_i {7};
  static constexpr unsigned int m_j {4};
  global_hits_alloc = m_i * longestamplicon + m_j + 1;
  global_hits_data = static_cast<unsigned int *>
    (xmalloc(global_hits_alloc * sizeof(unsigned int)));

  /* compute hash for all amplicons and store them in a hash table */

  const uint64_t hashtablesize {hash_alloc(amplicons)};
  bloom_a = bloom_init(hashtablesize);

  duplicates_found = 0;

  progress_init("Hashing sequences:", amplicons);
  for(auto k = 0U; k < amplicons; k++)
    {
      struct ampinfo_s * bp = ampinfo + k;
      bp->generation = 0;
      bp->swarmid = no_swarm;
      bp->next = no_swarm;
      bp->graft_cand = no_swarm;
      hash_insert(k);
      progress_update(k);
      if (duplicates_found != 0U) {
        break;
      }
    }

  if (duplicates_found != 0U)
    {
      fatal(error_prefix,
            "some fasta entries have identical sequences.\n"
            "Swarm expects dereplicated fasta files.\n"
            "Such files can be produced with swarm or vsearch:\n"
            " swarm -d 0 -w derep.fasta -o /dev/null input.fasta\n"
            "or\n"
            " vsearch --derep_fulllength input.fasta --sizein --sizeout --output derep.fasta\n");
    }

  progress_done();

  unsigned char * dir = nullptr;
  uint64_t * hearray = nullptr;

  /* for all amplicons, generate list of matching amplicons */

  network = static_cast<unsigned int*>
    (xmalloc(network_alloc * sizeof(unsigned int)));

  network_count = 0;

  pthread_mutex_init(&network_mutex, nullptr);
  network_amp = 0;
  progress_init("Building network: ", amplicons);
  auto * network_tr = new ThreadRunner(static_cast<int>(opt_threads),
                                       network_thread);
  network_tr->run();
  delete network_tr;
  pthread_mutex_destroy(&network_mutex);

  progress_done();


  /* dump network to file */
  if (not parameters.opt_network_file.empty()) {
    write_network_file(network_count, parameters);
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
          std::qsort(global_hits_data, global_hits_count,
                sizeof(unsigned int), compare_amp);

          /* add subseeds on list to current swarm */
          for(auto i = 0U; i < global_hits_count; i++) {
            add_amp_to_swarm(global_hits_data[i]);
          }

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
              std::qsort(global_hits_data, global_hits_count,
                    sizeof(unsigned int), compare_amp);

              /* add them to the swarm */
              for(auto i = 0U; i < global_hits_count; i++) {
                add_amp_to_swarm(global_hits_data[i]);
              }

              /* start with most abundant amplicon of next generation */
              if (global_hits_count != 0U) {
                subseed = global_hits_data[0];
              }
              else {
                subseed = no_swarm;
              }
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
          if (swarmsize > largest) {
            largest = swarmsize;
          }
          if (swarm_maxgen > maxgen) {
            maxgen = swarm_maxgen;
          }

          ++swarmcount;
        }
      progress_update(seed+1);
    }
  progress_done();

  if (global_hits_data != nullptr) {
    xfree(global_hits_data);
  }
  global_hits_data = nullptr;

  if (network != nullptr) {
    xfree(network);
  }
  network = nullptr;

  swarmcount_adjusted = swarmcount;

  /* fastidious */

  if (parameters.opt_fastidious)
    {
      std::fprintf(logfile, "\n");
      std::fprintf(logfile, "Results before fastidious processing:\n");
      std::fprintf(logfile, "Number of swarms:  %u\n", swarmcount);
      std::fprintf(logfile, "Largest swarm:     %u\n", largest);
      std::fprintf(logfile, "\n");

      uint64_t small_clusters = 0;
      uint64_t amplicons_in_small_clusters = 0;
      uint64_t nucleotides_in_small_clusters = 0;

      progress_init("Counting amplicons in heavy and light swarms",
                    swarmcount);

      for(auto i = 0ULL; i < swarmcount; i++)
        {
          struct swarminfo_s * sp = swarminfo + i;
          if (sp->mass < static_cast<uint64_t>(opt_boundary))
            {
              amplicons_in_small_clusters += sp->size;
              nucleotides_in_small_clusters += sp->sumlen;
              ++small_clusters;
            }
          progress_update(i+1);
        }
      progress_done();

      const uint64_t amplicons_in_large_clusters = amplicons - amplicons_in_small_clusters;
      const uint64_t large_clusters = swarmcount - small_clusters;

      std::fprintf(logfile, "Heavy swarms: %" PRIu64 ", with %" PRIu64 " amplicons\n",
              large_clusters, amplicons_in_large_clusters);
      std::fprintf(logfile, "Light swarms: %" PRIu64 ", with %" PRIu64 " amplicons\n",
              small_clusters, amplicons_in_small_clusters);
      std::fprintf(logfile, "Total length of amplicons in light swarms: %" PRIu64 "\n",
              nucleotides_in_small_clusters);

      if ((small_clusters == 0) or (large_clusters == 0))
        {
          std::fprintf(logfile, "Only light or heavy swarms found - "
                  "no need for further analysis.\n");
        }
      else
        {
          /* m: total size of Bloom filter in bits */
          /* k: number of hash functions */
          /* n: number of entries in the bloom filter */
          /* here: k=11 and m/n=18, that is 16 bits/entry */

          static constexpr unsigned int microvariants {7};
          static constexpr double hash_functions_per_bit {4.0 / 10};
          assert(parameters.opt_bloom_bits <= 64);  // larger than expected"
          assert(parameters.opt_bloom_bits >= 2);  // smaller than expected"
          auto bits = static_cast<unsigned int>(parameters.opt_bloom_bits);  // safe if opt_bloom_bits < UINT_MAX

          // int64_t k = int(bits * 0.693);    /* 11 */
          auto k =
            static_cast<unsigned int>(hash_functions_per_bit * bits); /* 6 */
          if (k < 1) {
            k = 1;
          }

          uint64_t m = bits * microvariants * nucleotides_in_small_clusters;
          static constexpr unsigned int min_total_bloom_filter_length_in_bits {64};

          const uint64_t memtotal = arch_get_memtotal();
          const uint64_t memused = arch_get_memused();

          if (parameters.opt_ceiling != 0)
            {
              const uint64_t memrest
                = one_megabyte * static_cast<uint64_t>(parameters.opt_ceiling) - memused;
              const auto new_bits =
                static_cast<unsigned int>(sizeof(uint64_t) * memrest / (microvariants * nucleotides_in_small_clusters));
              if (new_bits < bits)
                {
                  if (new_bits < 2) {
                    fatal(error_prefix, "Insufficient memory remaining for Bloom filter.");
                  }
                  std::fprintf(logfile, "Reducing memory used for Bloom filter due to --ceiling option.\n");
                  bits = new_bits;
                  // k = int(bits * 0.693);
                  k = static_cast<unsigned int>(hash_functions_per_bit * bits);
                  if (k < 1) {
                    k = 1;
                  }

                  m = bits * microvariants * nucleotides_in_small_clusters;
                }
            }

          if (m < min_total_bloom_filter_length_in_bits) {  // refactor C++17: std::clamp()
            m = min_total_bloom_filter_length_in_bits;  // at least 64 bits
          }

          if (memused + m / sizeof(uint64_t) > memtotal)
            {
              std::fprintf(logfile, "WARNING: Memory usage will probably exceed total amount of memory available.\n");
              std::fprintf(logfile, "Try to reduce memory footprint using the --bloom-bits or --ceiling options.\n");
            }

          std::fprintf(logfile,
                  "Bloom filter: bits=%u, m=%" PRIu64 ", k=%u, size=%.1fMB\n",
                  bits, m, k, static_cast<double>(m) / (sizeof(uint64_t) * one_megabyte));


          // m is in bits (divide by 8 to get bytes)
          // m is guaranteed to be at least 64 (see code above)
          static constexpr auto n_bits_in_a_byte = 8U;
          assert(m != 0);  // safeguard for future changes
          assert(m >= 64);
          const uint64_t n_bytes = ((m - 1) / n_bits_in_a_byte) + 1;
          bloom_f = bloomflex_init(n_bytes, k);


          /* Empty the old hash and bloom filter
             before we reinsert only the light swarm amplicons */

          hash_zap(hashtablesize);
          bloom_zap(bloom_a);

          progress_init("Adding light swarm amplicons to Bloom filter",
                        amplicons_in_small_clusters);

          /* process amplicons in order from least to most abundant */
          /* but stop when all amplicons in small clusters are processed */

          light_variants = 0;

          pthread_mutex_init(&light_mutex, nullptr);
          light_progress = 0;
          light_amplicon_count = amplicons_in_small_clusters;
          light_amplicon = amplicons - 1;
          auto * tr = new ThreadRunner(static_cast<int>(opt_threads),
                                       mark_light_thread);
          tr->run();
          delete tr;
          pthread_mutex_destroy(&light_mutex);

          progress_done();

          std::fprintf(logfile,
                  "Generated %" PRIu64 " variants from light swarms\n",
                  light_variants);

          progress_init("Checking heavy swarm amplicons against Bloom filter",
                        amplicons_in_large_clusters);

          /* process amplicons in order from most to least abundant */
          /* but stop when all amplicons in large clusters are processed */

          pthread_mutex_init(&graft_mutex, nullptr);

          heavy_variants = 0;

          pthread_mutex_init(&heavy_mutex, nullptr);
          heavy_progress = 0;
          heavy_amplicon_count = amplicons_in_large_clusters;
          heavy_amplicon = 0;
          auto * heavy_tr
            = new ThreadRunner(static_cast<int>(opt_threads),
                               check_heavy_thread);
          heavy_tr->run();
          delete heavy_tr;
          pthread_mutex_destroy(&heavy_mutex);

          progress_done();

          bloomflex_exit(bloom_f);

          pthread_mutex_destroy(&graft_mutex);

          std::fprintf(logfile, "Heavy variants: %" PRIu64 "\n", heavy_variants);
          std::fprintf(logfile, "Got %" PRId64 " graft candidates\n", graft_candidates);
          const unsigned int grafts = attach_candidates(amplicons);
          std::fprintf(logfile, "Made %u grafts\n", grafts);
          std::fprintf(logfile, "\n");
        }
    }


  /* dump swarms */
  if (parameters.opt_mothur) {
    write_swarms_mothur_format(swarmcount, parameters);
  }
  else {
    write_swarms_default_format(swarmcount, parameters);
  }

  /* dump seeds in fasta format with sum of abundances */
  if (not parameters.opt_seeds.empty()) {
    write_representative_sequences(swarmcount, parameters);
  }

  /* output internal structure */
  if (not parameters.opt_internal_structure.empty()) {
    write_structure_file(swarmcount, parameters);
  }

  /* output swarm in uclust format */
  if (uclustfile != nullptr) {
    write_swarms_uclust_format(swarmcount, parameters, dir, hearray);
  }

  /* output statistics to file */
  if (statsfile != nullptr) {
    write_stats_file(swarmcount, parameters);
  }

  std::fprintf(logfile, "\n");
  std::fprintf(logfile, "Number of swarms:  %" PRIu64 "\n", swarmcount_adjusted);
  std::fprintf(logfile, "Largest swarm:     %u\n", largest);
  std::fprintf(logfile, "Max generations:   %u\n", maxgen);

  bloom_exit(bloom_a);
  hash_free();

  if (swarminfo != nullptr) {
    xfree(swarminfo);
  }

  delete [] ampinfo;
  ampinfo = nullptr;

#ifdef HASHSTATS
  std::fprintf(logfile, "Tries:      %12lu\n", tries);
  std::fprintf(logfile, "Bloom m:    %12lu\n", bloom_matches);
  std::fprintf(logfile, "Hits:       %12lu\n", hits);
  std::fprintf(logfile, "Success:    %12lu\n", success);
  std::fprintf(logfile, "Bingo:      %12lu\n", bingo);
  std::fprintf(logfile, "Collisions: %12lu\n", collisions);
#endif
}
