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
#include <algorithm>  // std::sort()
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <climits>  // UINT_MAX
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc(), size_t
#include <cstdlib>  // qsort()
#include <cstring>  // std::memcmp
#include <limits>  // unsigned int max
#include <memory>  // unique pointer
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
constexpr unsigned int no_swarm {std::numeric_limits<unsigned int>::max()};

static uint64_t duplicates_found {0};  // several function calls

/* Information about each amplicon */

struct ampinfo_s
{
  unsigned int swarmid {no_swarm};
  unsigned int parent {0U};
  unsigned int generation {0U};
  unsigned int next {no_swarm};        /* amp id of next amplicon in swarm */
  unsigned int graft_cand {no_swarm};  /* amp id of potential grafting parent (fastid.) */
  unsigned int link_start {0U};
  unsigned int link_count {0U};
};

// refactoring: can't be eliminated yet, because of the pthread barrier
static struct ampinfo_s * ampinfo = nullptr;

/* Information about each swarm (cluster) */

struct swarminfo_s
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
};  // total of 40 bytes (five 64-bit machine words)

static struct swarminfo_s* swarminfo = nullptr;

struct graft_cand
{
  unsigned int parent;
  unsigned int child;
};

/* Information about potential grafts */
static int64_t graft_candidates {0};
static pthread_mutex_t graft_mutex;

static unsigned int current_swarm_tail {0};

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


inline auto check_amp_identical(unsigned int amp1,
                                unsigned int amp2) -> bool
{
  /* amplicon are identical if they have the same length, and the
     exact same sequence */
  const auto amp1_seqlen = db_getsequencelen(amp1);
  const auto amp2_seqlen = db_getsequencelen(amp2);

  return ((amp1_seqlen == amp2_seqlen) and
          (std::memcmp(db_getsequence(amp1),
                  db_getsequence(amp2),
                  nt_bytelength(amp1_seqlen)) == 0));
}


inline auto hash_insert(unsigned int amp) -> void
{
  /* find the first empty bucket */
  const auto hash = db_gethash(amp);
  auto index = hash_getindex(hash);
  auto duplicate = false;
  while (hash_is_occupied(index))
    {
      if (hash_compare_value(index, hash) and
          check_amp_identical(amp, hash_get_data(index))) {
        duplicate = true;
      }
      index = hash_getnextindex(index);
    }

  if (duplicate) {
    ++duplicates_found;
  }

  hash_set_occupied(index);
  hash_set_value(index, hash);
  hash_set_data(index, amp);

  bloom_set(bloom_a, hash);
}


/******************** FASTIDIOUS START ********************/


auto attach(unsigned int seed, unsigned int amp,
            std::vector<struct ampinfo_s> & ampinfo_v,
            std::vector<struct swarminfo_s> & swarminfo_v) -> void
{
  /* graft light swarm (amp) on heavy swarm (seed) */

  swarminfo_s & heavy_swarm = swarminfo_v[ampinfo_v[seed].swarmid];
  swarminfo_s & light_swarm = swarminfo_v[ampinfo_v[amp].swarmid];

  // attach the seed of the light swarm to the tail of the heavy swarm (refactoring: unclear)
  ampinfo_v[heavy_swarm.last].next = light_swarm.seed;
  heavy_swarm.last = light_swarm.last;

  // Update swarm info
  heavy_swarm.size += light_swarm.size;
  heavy_swarm.singletons += light_swarm.singletons;
  heavy_swarm.mass += light_swarm.mass;
  heavy_swarm.sumlen += light_swarm.sumlen;
  /* maxgen is untouched */

  /* flag attachment to avoid doing it again */
  light_swarm.attached = true;

  // Update overall stats
  if (heavy_swarm.size > largest) {
    largest = heavy_swarm.size;
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


// refactoring: replace with std::count_if()
auto count_pairs(unsigned int const amplicon_count,
                 std::vector<struct ampinfo_s> & ampinfo_v) -> unsigned int {
  auto counter = 0U;
  for(auto i = 0U; i < amplicon_count; i++) {
    if (ampinfo_v[i].graft_cand != no_swarm) {
      ++counter;
    }
  }
  return counter;
}


auto attach_candidates(unsigned int amplicon_count,
                       std::vector<struct ampinfo_s> & ampinfo_v,
                       std::vector<struct swarminfo_s> & swarminfo_v) -> unsigned int
{
  auto const pair_count = count_pairs(amplicon_count, ampinfo_v);

  progress_init("Grafting light swarms on heavy swarms", pair_count);

  /* allocate memory */
  std::vector<struct graft_cand> graft_array(pair_count);

  /* fill in */
  auto ticker = 0U;  // refactoring: replace with a transform algorithm
  for(auto i = 0U; i < amplicon_count; i++) {
    if (ampinfo_v[i].graft_cand != no_swarm)
      {
        graft_array[ticker].parent = ampinfo_v[i].graft_cand;
        graft_array[ticker].child = i;  // so two children cannot have the same uint value
        ++ticker;
      }
  }

  /* sort */
  auto compare_grafts = [](struct graft_cand const& lhs,
                           struct graft_cand const& rhs) -> bool {
    // sort by parent index (lowest index first)
    if (lhs.parent < rhs.parent) {
      return true;
    }
    if (lhs.parent > rhs.parent) {
      return false;
    }
    // ...then ties are sorted by child index (lowest index first)
    // assert(lhs.child != rhs.child); // refactoring: should be true by construction! It is not, why?
    // assert(lhs.child < rhs.child); // refactoring: same as above!
    if (lhs.child < rhs.child) {
      return true;
    }
    return false;
  };

  std::sort(graft_array.begin(), graft_array.end(), compare_grafts);

  /* attach in order */
  auto grafts = 0U;
  for(auto i = 0U; i < pair_count; i++)
    {
      const auto parent = graft_array[i].parent;
      const auto child  = graft_array[i].child;

      if (swarminfo_v[ampinfo_v[child].swarmid].attached)
        {
          /* this light swarm is already attached */
          ampinfo_v[child].graft_cand = no_swarm;
        }
      else
        {
          /* attach child to parent */
          attach(parent, child, ampinfo_v, swarminfo_v);
          ++grafts;
        }
      progress_update(i + 1);
    }
  progress_done();
  return grafts;
}


auto hash_check_attach(char * seed_sequence,
                       unsigned int seed_seqlen,
                       struct var_s & var,
                       unsigned int seed) -> bool
{
  /* seed is the original large swarm seed */

  /* compute hash and corresponding hash table index */
  const auto hash = var.hash;
  auto index = hash_getindex(hash);

  /* find matching buckets */

  while (hash_is_occupied(index))
    {
      if (hash_compare_value(index, hash))
        {
          /* check that mass is below threshold */
          const auto amp = hash_get_data(index);

          /* make absolutely sure sequences are identical */
          auto *amp_sequence = db_getsequence(amp);
          const auto amp_seqlen = db_getsequencelen(amp);
          if (check_variant(seed_sequence, seed_seqlen, var, amp_sequence, amp_seqlen))
            {
              add_graft_candidate(seed, amp);
              return true;
            }
        }
      index = hash_getnextindex(index);
    }
  return false;
}


inline auto check_heavy_var_2(std::vector<char>& seq,
                              unsigned int seqlen,
                              unsigned int seed,
                              std::vector<struct var_s>& variant_list) -> uint64_t
{
  /* Check second generation microvariants of the heavy swarm amplicons
     and see if any of them are identical to a light swarm amplicon. */

  uint64_t matches = 0;

  const auto hash = zobrist_hash(reinterpret_cast<unsigned char *>(seq.data()), seqlen);
  const auto variant_count = generate_variants(seq.data(), seqlen, hash, variant_list);  // refactoring: seq.data() not fixable while db returns char*

  for(auto i = 0U; i < variant_count; i++) {
    if (bloom_get(bloom_a, variant_list[i].hash) and
        hash_check_attach(seq.data(), seqlen, variant_list[i], seed)) {
      ++matches;
    }
  }

  return matches;
}


auto check_heavy_var(struct bloomflex_s * bloom,
                     std::vector<char>& varseq,
                     unsigned int seed,
                     uint64_t & number_of_matches,
                     uint64_t & number_of_variants,
                     std::vector<struct var_s>& variant_list,
                     std::vector<struct var_s>& variant_list2) -> void
{
  /*
    bloom is a bloom filter in which to check the variants
    varseq is a buffer large enough to hold any sequence + 1 insertion
    seed is the original seed
    number_of_matches is where to store number of matches
    number_of_variants is where to store number of variants
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

  uint64_t matches = 0;

  auto *sequence = db_getsequence(seed);
  const auto seqlen = db_getsequencelen(seed);
  const auto hash = db_gethash(seed);
  const auto variant_count = generate_variants(sequence, seqlen, hash, variant_list);

  for(auto i = 0U; i < variant_count; i++)
    {
      struct var_s & var = variant_list[i];
      if (bloomflex_get(bloom, var.hash))
        {
          auto varlen = 0U;
          generate_variant_sequence(sequence, seqlen,
                                    var, varseq, varlen);
          matches += check_heavy_var_2(varseq,
                                       varlen,
                                       seed,
                                       variant_list2);
        }
    }

  number_of_matches = matches;
  number_of_variants = variant_count;
}


auto check_heavy_thread(int64_t t) -> void
{
  static constexpr auto i = 7U;  // max number of microvariants = 7 * len + 4
  static constexpr auto j = 4U;  //                               i * len + j
  static constexpr auto nt_per_uint64 = 32U;  // 32 nucleotides can fit in a uint64
  (void) t;

  std::vector<struct var_s> variant_list(i * longestamplicon + j);
  std::vector<struct var_s> variant_list2(i * (longestamplicon + 1) + j);

  const std::size_t size =
    sizeof(uint64_t) * ((db_getlongestsequence() + 2 + nt_per_uint64 - 1) / nt_per_uint64);
  std::vector<char> buffer1(size);
  pthread_mutex_lock(&heavy_mutex);
  while ((heavy_amplicon < amplicons) and
         (heavy_progress < heavy_amplicon_count))
    {
      const auto heavy_amplicon_id = heavy_amplicon++;
      if (swarminfo[ampinfo[heavy_amplicon_id].swarmid].mass >=
          static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++heavy_progress);
          pthread_mutex_unlock(&heavy_mutex);
          uint64_t number_of_matches {0};
          uint64_t number_of_variants {0};
          check_heavy_var(bloom_f, buffer1, heavy_amplicon_id,
                          number_of_matches, number_of_variants,
                          variant_list, variant_list2);
          pthread_mutex_lock(&heavy_mutex);
          heavy_variants += number_of_variants;
        }
    }
  pthread_mutex_unlock(&heavy_mutex);
}


auto mark_light_var(struct bloomflex_s * bloom,
                        unsigned int seed,
                        std::vector<struct var_s>& variant_list) -> uint64_t
{
  /*
    add all microvariants of seed to Bloom filter

    bloom is a BloomFilter in which to enter the variants
    seed is the original seed
  */

  hash_insert(seed);

  auto *sequence = db_getsequence(seed);
  const auto seqlen = db_getsequencelen(seed);
  const auto hash = db_gethash(seed);
  const auto variant_count = generate_variants(sequence, seqlen, hash, variant_list);

  for(auto i = 0U; i < variant_count; i++) {
    bloomflex_set(bloom, variant_list[i].hash);
  }

  return variant_count;
}


void mark_light_thread(int64_t t)
{
  static constexpr auto i = 7U;  // max number of microvariants = 7 * len + 4
  static constexpr auto j = 4U;  //                               i * len + j

  (void) t;

  std::vector<struct var_s> variant_list(i * longestamplicon + j);

  pthread_mutex_lock(&light_mutex);
  while (light_progress < light_amplicon_count)
    {
      const auto light_amplicon_id = light_amplicon--;
      if (swarminfo[ampinfo[light_amplicon_id].swarmid].mass <
          static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++light_progress);
          pthread_mutex_unlock(&light_mutex);
          const auto variant_count = mark_light_var(bloom_f, light_amplicon_id,
                                                        variant_list);
          pthread_mutex_lock(&light_mutex);
          light_variants += variant_count;
        }
    }
  pthread_mutex_unlock(&light_mutex);
}


/******************** FASTIDIOUS END ********************/


inline auto find_variant_matches(unsigned int seed,
                                 struct var_s & var,
                                 std::vector<unsigned int>& hits_data,
                                 unsigned int & hits_count) -> void
{
  if (not bloom_get(bloom_a, var.hash)) {
    return;
  }

  /* compute hash and corresponding hash table index */

  auto index = hash_getindex(var.hash);

  /* find matching buckets */

  while (hash_is_occupied(index))
    {
      if (hash_compare_value(index, var.hash))
        {
          const auto amp = hash_get_data(index);

          /* avoid self */
          if (seed != amp) {
            if ((opt_no_otu_breaking) or
                (db_getabundance(seed) >= db_getabundance(amp)))
              {
                auto *seed_sequence = db_getsequence(seed);
                const auto seed_seqlen = db_getsequencelen(seed);

                auto *amp_sequence = db_getsequence(amp);
                const auto amp_seqlen = db_getsequencelen(amp);

                if (check_variant(seed_sequence, seed_seqlen,
                                  var,
                                  amp_sequence, amp_seqlen))
                  {
                    hits_data[hits_count++] = amp;
                    break;
                  }
              }
          }
        }
      index = hash_getnextindex(index);
    }
}


auto check_variants(unsigned int seed,
                    std::vector<struct var_s> & variant_list,
                    std::vector<unsigned int>& hits_data) -> unsigned int
{
  auto hits_count = 0U;

  auto *sequence = db_getsequence(seed);
  const auto seqlen = db_getsequencelen(seed);
  const auto hash = db_gethash(seed);
  const auto variant_count = generate_variants(sequence, seqlen, hash, variant_list);

  // refactoring: range-based for loop over variant_list truncated to variant_count?
  for(auto i = 0U; i < variant_count; i++) {
    find_variant_matches(seed, variant_list[i], hits_data, hits_count);
  }

  return hits_count;
}


auto network_thread(int64_t t) -> void
{
  static constexpr auto i = 7U;  // max number of microvariants = 7 * len + 4
  static constexpr auto j = 4U;  //                               i * len + j

  (void) t;

  std::vector<unsigned int> hits_data(i * longestamplicon + j + 1);
  std::vector<struct var_s> variant_list(i * longestamplicon + j + 1);

  pthread_mutex_lock(&network_mutex);
  while (network_amp < amplicons)
    {
      const auto amp = network_amp++;
      progress_update(amp);

      pthread_mutex_unlock(&network_mutex);

      const auto hits_count = check_variants(amp, variant_list, hits_data);
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
}


auto process_seed(unsigned int seed,
                  std::vector<struct ampinfo_s> & ampinfo_v,
                  std::vector<unsigned int> & global_hits_v,
                  unsigned int & global_hits_count) -> void
{
  /* update swarm stats */
  struct ampinfo_s const & bp = ampinfo_v[seed];

  ++swarmsize;
  if (bp.generation > swarm_maxgen) {
    swarm_maxgen = bp.generation;
  }
  const auto abundance = db_getabundance(seed);
  abundance_sum += abundance;
  if (abundance == 1) {
    ++singletons;
  }
  swarm_sumlen += db_getsequencelen(seed);

  const auto link_start = ampinfo_v[seed].link_start;
  const auto link_count = ampinfo_v[seed].link_count;
  auto global_hits_alloc = global_hits_v.size();

  if (global_hits_count + link_count > global_hits_alloc)
    {
      while (global_hits_count + link_count > global_hits_alloc) {
        global_hits_alloc += 4 * one_kilobyte;
      }
      global_hits_v.resize(global_hits_alloc);
      global_hits_data = global_hits_v.data();
    }

  for(auto offset = 0U; offset < link_count; offset++)
    {
      const auto amp = network[link_start + offset];

      if (ampinfo_v[amp].swarmid == no_swarm)
        {
          global_hits_v[global_hits_count] = amp;
          ++global_hits_count;

          /* update info */
          ampinfo_v[amp].swarmid = ampinfo_v[seed].swarmid;
          ampinfo_v[amp].generation = ampinfo_v[seed].generation + 1;
          ampinfo_v[amp].parent = seed;
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
  auto status = 1;  // default is *x > *y

  assert(*x != *y);  // '*x == *y' is not expected at that stage

  if (*x < *y) {
    status = -1;
  }

  return status;
}


inline auto add_amp_to_swarm(unsigned int const amp,
                             std::vector<struct ampinfo_s> & ampinfo_v) -> void
{
  /* add to swarm */
  ampinfo_v[current_swarm_tail].next = amp;
  current_swarm_tail = amp;
}


auto write_network_file(const unsigned int number_of_networks,
                        struct Parameters const & parameters,
                        std::vector<struct ampinfo_s> & ampinfo_v) -> void {
  // network = cluster with at least two sequences (no singletons)
  progress_init("Dumping network:  ", number_of_networks);

  uint64_t n_processed = 0;  // refactoring: reduce scope (add to for loop init)
  assert(ampinfo_v.size() == amplicons);
  auto counter = 0ULL;
  for(auto const& amplicon: ampinfo_v) {
      const auto link_start = amplicon.link_start;
      const auto link_count = amplicon.link_count;

      std::qsort(network + link_start,
            link_count,
            sizeof(unsigned int),
            compare_amp);

      for(auto link = 0U; link < link_count; link++)
        {
          const auto neighbour = network[link_start + link];
          fprint_id(network_file, counter, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(network_file, "\t");
          fprint_id(network_file, neighbour, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(network_file, "\n");
          ++n_processed;
        }
      progress_update(n_processed);
      ++counter;
    }
  progress_done();
}


auto write_swarms_default_format(const unsigned int swarmcount,
                                 struct Parameters const & parameters,
                                 std::vector<struct ampinfo_s> & ampinfo_v,
                                 std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  progress_init("Writing swarms:   ", swarmcount);

  for(auto i = 0U; i < swarmcount; i++) {
    if (swarminfo_v[i].attached) {
      continue;
    }

    const auto seed = swarminfo_v[i].seed;
    for(auto a = seed; a != no_swarm; a = ampinfo_v[a].next) {
      if (a != seed) {
        std::fputc(sepchar, outfile);
      }
      fprint_id(outfile, a,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    }
    std::fputc('\n', outfile);
    progress_update(i + 1);
  }

  progress_done();
}


auto write_swarms_mothur_format(const unsigned int swarmcount,
                                struct Parameters const & parameters,
                                std::vector<struct ampinfo_s> & ampinfo_v,
                                std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  progress_init("Writing swarms:   ", swarmcount);

  std::fprintf(outfile, "swarm_%" PRId64 "\t%" PRIu64,
          parameters.opt_differences, swarmcount_adjusted);

  for(auto i = 0U; i < swarmcount; i++) {
    if (swarminfo_v[i].attached) {
      continue;
    }

    const auto seed = swarminfo_v[i].seed;
    for(auto a = seed; a != no_swarm; a = ampinfo_v[a].next) {
      if (a == seed) {
        std::fputc('\t', outfile);
      }
      else {
        std::fputc(',', outfile);
      }
      fprint_id(outfile, a,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    }
    progress_update(i + 1);
  }

  std::fputc('\n', outfile);

  progress_done();
}


auto write_swarms_uclust_format(const unsigned int swarmcount,
                                struct Parameters const & parameters,
                                std::vector<struct ampinfo_s> & ampinfo_v,
                                std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  static constexpr auto one_hundred = 100U;
  auto cluster_no = 0U;
  const auto score_matrix_63 = create_score_matrix<int64_t>(parameters.penalty_mismatch);
  std::vector<unsigned char> directions(longestamplicon * longestamplicon);
  std::vector<uint64_t> hearray(2 * longestamplicon);

  progress_init("Writing UCLUST:   ", swarmcount);

  for(auto swarmid = 0U; swarmid < swarmcount; swarmid++)
    {
      if (swarminfo_v[swarmid].attached) {
        continue;
      }

      const auto seed = swarminfo_v[swarmid].seed;

      struct ampinfo_s const & bp = ampinfo_v[seed];

      std::fprintf(uclustfile, "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                   cluster_no,
                   swarminfo_v[swarmid].size);
      fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(uclustfile, "\t*\n");

      std::fprintf(uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                   cluster_no,
                   db_getsequencelen(seed));
      fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(uclustfile, "\t*\n");

      for(auto a = bp.next; a != no_swarm; a = ampinfo_v[a].next)
        {
          auto * dseq = db_getsequence(a);
          const auto dlen = db_getsequencelen(a);
          auto * qseq = db_getsequence(seed);
          const auto qlen = db_getsequencelen(seed);

          uint64_t nwdiff = 0;
          char * nwalignment = nullptr;  // CIGAR string
          uint64_t nwalignmentlength = 0;

          nw(dseq, dlen, qseq, qlen,
             score_matrix_63, static_cast<unsigned long int>(penalty_gapopen),
             static_cast<unsigned long int>(penalty_gapextend),
             nwdiff, nwalignmentlength, nwalignment,
             directions, hearray);

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
      progress_update(swarmid);
    }
  progress_done();
}


auto write_representative_sequences(const unsigned int swarmcount,
                                    struct Parameters const & parameters,
                                    std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  progress_init("Writing seeds:    ", swarmcount);

  std::vector<unsigned int> sorter(swarmcount);
  for(auto i = 0U; i < swarmcount; i++) {
    sorter[i] = i;
  }

  auto compare_mass = [&swarminfo_v](unsigned int const lhs,
                                     unsigned int const rhs) -> bool
  {
    const swarminfo_s & swarm_x = swarminfo_v[lhs];
    const swarminfo_s & swarm_y = swarminfo_v[rhs];

    const auto mass_x = swarm_x.mass;
    const auto mass_y = swarm_y.mass;

    // sort seeds by decreasing mass
    if (mass_x > mass_y) {
      return true;
    }
    if (mass_x < mass_y) {
      return false;
    }
    // ...then ties are sorted by headers (alphabetical order)
    const auto status = std::strcmp(db_getheader(swarm_x.seed),
                                    db_getheader(swarm_y.seed));
    if (status < 0) {
      return true;
    }
    // assert(status != 0); // all headers are unique
    return false;
  };

  std::sort(sorter.begin(), sorter.end(), compare_mass);

  for(auto j = 0U; j < swarmcount; j++)
    {
      const auto index = sorter[j];
      if (swarminfo_v[index].attached) {
        continue;
      }
      const auto seed = swarminfo_v[index].seed;
      std::fprintf(fp_seeds, ">");
      fprint_id_with_new_abundance(fp_seeds, seed, swarminfo_v[index].mass, parameters.opt_usearch_abundance);
      std::fprintf(fp_seeds, "\n");
      db_fprintseq(fp_seeds, seed);
      progress_update(index + 1);
    }

  progress_done();
}


auto write_structure_file(const unsigned int swarmcount,
                          struct Parameters const & parameters,
                          std::vector<struct ampinfo_s> & ampinfo_v,
                          std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  auto cluster_no = 0U;

  progress_init("Writing structure:", swarmcount);

  for(auto swarmid = 0U; swarmid < swarmcount; swarmid++)
    {
      if (swarminfo_v[swarmid].attached) {
        continue;
      }
      const auto seed = swarminfo_v[swarmid].seed;

      struct ampinfo_s const & bp = ampinfo_v[seed];

      for(auto a = bp.next; a != no_swarm; a = ampinfo_v[a].next)
        {
          const auto graft_parent = ampinfo_v[a].graft_cand;
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
                           ampinfo_v[graft_parent].generation + 1);
            }

          const auto parent = ampinfo_v[a].parent;
          if (parent != no_swarm)
            {
              fprint_id_noabundance(internal_structure_file, parent, parameters.opt_usearch_abundance);
              std::fprintf(internal_structure_file, "\t");
              fprint_id_noabundance(internal_structure_file, a, parameters.opt_usearch_abundance);
              std::fprintf(internal_structure_file,
                           "\t%u\t%u\t%u\n",
                           1U,
                           cluster_no + 1,
                           ampinfo_v[a].generation);
            }
        }

      ++cluster_no;
      progress_update(swarmid);
    }
  progress_done();
}


auto write_stats_file(const unsigned int swarmcount,
                      struct Parameters const & parameters,
                      std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  progress_init("Writing stats:    ", swarmcount);
  for(auto i = 0ULL; i < swarmcount; i++)
    {
      struct swarminfo_s const & sp = swarminfo_v[i];
      if (sp.attached) {
        continue;
      }
      std::fprintf(statsfile, "%u\t%" PRIu64 "\t", sp.size, sp.mass);
      fprint_id_noabundance(statsfile, sp.seed, parameters.opt_usearch_abundance);
      std::fprintf(statsfile, "\t%" PRIu64 "\t%u\t%u\t%u\n",
                   db_getabundance(sp.seed),
                   sp.singletons, sp.maxgen, sp.maxgen);
      progress_update(i);
    }
  progress_done();
}


auto algo_d1_run(struct Parameters const & parameters) -> void
{
  longestamplicon = db_getlongestsequence();
  amplicons = db_getsequencecount();

  std::vector<struct ampinfo_s> ampinfo_v(amplicons);
  ampinfo = ampinfo_v.data();

  std::vector<struct swarminfo_s> swarminfo_v(one_kilobyte);
  swarminfo = swarminfo_v.data();

  // max number of microvariants = 7 * len + 4
  static constexpr auto m_i = 7U;
  static constexpr auto m_j = 4U;
  const auto global_hits_alloc = m_i * longestamplicon + m_j + 1;
  std::vector<unsigned int> global_hits_v(global_hits_alloc);
  global_hits_data = global_hits_v.data();

  /* compute hash for all amplicons and store them in a hash table */

  const uint64_t hashtablesize {hash_alloc(amplicons)};
  bloom_a = bloom_init(hashtablesize);

  duplicates_found = 0;

  progress_init("Hashing sequences:", amplicons);
  for(auto k = 0U; k < amplicons; k++)
    {
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

  /* for all amplicons, generate list of matching amplicons */

  // refactoring: not possible because of the pthread barrier
  // std::vector<unsigned int> network_v(network_alloc);
  // network = network_v.data();
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
    write_network_file(network_count, parameters, ampinfo_v);
  }


  /* for each non-swarmed amplicon look for subseeds ... */

  auto swarmcount = 0U;
  progress_init("Clustering:       ", amplicons);

  for(auto seed = 0U; seed < amplicons; seed++)
    {
      struct ampinfo_s & ap = ampinfo_v[seed];

      if (ap.swarmid == no_swarm)
        {
          /* start a new swarm with a new initial seed */

          ap.swarmid = swarmcount;
          ap.generation = 0;
          ap.parent = no_swarm;
          ap.next = no_swarm;

          /* link up this initial seed in the list of swarms */
          current_swarm_tail = seed;

          /* initialize swarm stats */
          swarmsize = 0;
          swarm_maxgen = 0;
          abundance_sum = 0;
          singletons = 0;
          swarm_sumlen = 0;

          /* init list */
          auto global_hits_count = 0U;

          /* find the first generation matches */
          process_seed(seed, ampinfo_v, global_hits_v, global_hits_count);

          /* sort hits */
          std::sort(global_hits_v.begin(), global_hits_v.begin() + global_hits_count);

          /* add subseeds on list to current swarm */
          for(auto i = 0U; i < global_hits_count; i++) {
            add_amp_to_swarm(global_hits_v[i], ampinfo_v);
          }

          /* find later generation matches */
          auto subseed = ap.next;
          while(subseed != no_swarm)
            {
              /* process all subseeds of this generation */
              global_hits_count = 0;

              while(subseed != no_swarm)
                {
                  process_seed(subseed, ampinfo_v, global_hits_v, global_hits_count);
                  subseed = ampinfo_v[subseed].next;
                }

              /* sort all of this generation */
              std::sort(global_hits_v.begin(), global_hits_v.begin() + global_hits_count);

              /* add them to the swarm */
              for(auto i = 0U; i < global_hits_count; i++) {
                add_amp_to_swarm(global_hits_v[i], ampinfo_v);
              }

              /* start with most abundant amplicon of next generation */
              if (global_hits_count != 0U) {
                subseed = global_hits_v[0];
              }
              else {
                subseed = no_swarm;
              }
            }

          if (swarmcount >= swarminfo_v.capacity())
            {
              /* allocate memory for more swarms... */
              // 1,024 times struct size, so 40,960 bytes
              swarminfo_v.resize(swarminfo_v.capacity() + one_kilobyte);
              swarminfo = swarminfo_v.data();
            }

          struct swarminfo_s & sp = swarminfo_v[swarmcount];

          sp.seed = seed;
          sp.size = swarmsize;
          sp.mass = abundance_sum;
          sp.sumlen = swarm_sumlen;
          sp.singletons = singletons;
          sp.maxgen = swarm_maxgen;
          sp.last = current_swarm_tail;
          sp.attached = false;

          /* update overall stats */
          if (swarmsize > largest) {
            largest = swarmsize;
          }
          if (swarm_maxgen > maxgen) {
            maxgen = swarm_maxgen;
          }

          ++swarmcount;
        }
      progress_update(seed + 1);
    }
  progress_done();

  // refactoring: trim vectors (remove allocated unused elements) (FAIL)
  // swarminfo_v.resize(swarmcount);
  // swarminfo_v.shrink_to_fit();

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

      // refactoring: move to function that returns a struct cluster_stats
      progress_init("Counting amplicons in heavy and light swarms",
                    swarmcount);

      for(auto i = 0ULL; i < swarmcount; i++)
        {
          struct swarminfo_s const & sp = swarminfo_v[i];
          if (sp.mass < static_cast<uint64_t>(opt_boundary))
            {
              amplicons_in_small_clusters += sp.size;
              nucleotides_in_small_clusters += sp.sumlen;
              ++small_clusters;
            }
          progress_update(i + 1);
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

          uint64_t m = nucleotides_in_small_clusters * microvariants * bits;
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

                  m = nucleotides_in_small_clusters * microvariants * bits;
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
          // refactoring: light_tr with unique ptr (see below for heavy_tr)
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
          {
            // refactoring C++14: use std::make_unique
            std::unique_ptr<ThreadRunner> heavy_tr (new ThreadRunner(static_cast<int>(opt_threads), check_heavy_thread));
            heavy_tr->run();
          }

          pthread_mutex_destroy(&heavy_mutex);

          progress_done();

          bloomflex_exit(bloom_f);

          pthread_mutex_destroy(&graft_mutex);

          std::fprintf(logfile, "Heavy variants: %" PRIu64 "\n", heavy_variants);
          std::fprintf(logfile, "Got %" PRId64 " graft candidates\n", graft_candidates);
          const unsigned int grafts = attach_candidates(amplicons, ampinfo_v, swarminfo_v);
          std::fprintf(logfile, "Made %u grafts\n", grafts);
          std::fprintf(logfile, "\n");
        }
    }


  /* dump swarms */
  if (parameters.opt_mothur) {
    write_swarms_mothur_format(swarmcount, parameters, ampinfo_v, swarminfo_v);
  }
  else {
    write_swarms_default_format(swarmcount, parameters, ampinfo_v, swarminfo_v);
  }

  /* dump seeds in fasta format with sum of abundances */
  if (not parameters.opt_seeds.empty()) {
    write_representative_sequences(swarmcount, parameters, swarminfo_v);
  }

  /* output internal structure */
  if (not parameters.opt_internal_structure.empty()) {
    write_structure_file(swarmcount, parameters, ampinfo_v, swarminfo_v);
  }

  /* output swarms in uclust format */
  if (uclustfile != nullptr) {
    write_swarms_uclust_format(swarmcount, parameters, ampinfo_v, swarminfo_v);
  }

  /* output statistics to file */
  if (statsfile != nullptr) {
    write_stats_file(swarmcount, parameters, swarminfo_v);
  }

  std::fprintf(logfile, "\n");
  std::fprintf(logfile, "Number of swarms:  %" PRIu64 "\n", swarmcount_adjusted);
  std::fprintf(logfile, "Largest swarm:     %u\n", largest);
  std::fprintf(logfile, "Max generations:   %u\n", maxgen);

  bloom_exit(bloom_a);
  hash_free();

  swarminfo = nullptr;
  ampinfo = nullptr;
}
