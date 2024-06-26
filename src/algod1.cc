/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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
#include "utils/cigar.h"
#include "utils/nt_codec.h"
#include "utils/opt_boundary.h"
#include "utils/opt_no_cluster_breaking.h"
#include "utils/progress.h"
#include "utils/score_matrix.h"
#include "utils/threads.h"
#include "zobrist.h"
#include <algorithm>  // std::sort(), std::reverse(), std::max()
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc(), size_t
#include <cstdlib>  // qsort()
#include <cstring>  // std::strcmp
#include <iterator>  // std::next
#include <limits>  // unsigned int max
#include <memory>  // unique pointer
#include <numeric>  // std::iota
#include <pthread.h>
#include <string>
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


constexpr unsigned int one_kilobyte {1U << 10U};  // 1,024 bytes
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

std::vector<unsigned int> network_v;
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
          std::equal(db_getsequence(amp1),
                     std::next(db_getsequence(amp1), nt_bytelength(amp1_seqlen)),
                     db_getsequence(amp2)));
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
  largest = std::max(heavy_swarm.size, largest);

  --swarmcount_adjusted;
}


auto add_graft_candidate(unsigned int seed, unsigned int amp) -> void
{
  pthread_mutex_lock(&graft_mutex);
  ++graft_candidates;
  assert(amp <= std::numeric_limits<std::ptrdiff_t>::max());
  auto const signed_position = static_cast<std::ptrdiff_t>(amp);
  auto & amplicon = *std::next(ampinfo, signed_position);
  // if there is no heavy candidate to graft amp, or if seed is
  // earlier in the sorting order, then we change the attachment to
  // seed
  if ((amplicon.graft_cand == no_swarm) or (amplicon.graft_cand > seed)) {
    amplicon.graft_cand = seed;
  }
  pthread_mutex_unlock(&graft_mutex);
}


// C++17 refactoring: replace with std::count_if()
auto count_pairs(unsigned int const amplicon_count,
                 std::vector<struct ampinfo_s> & ampinfo_v) -> unsigned int {
  auto counter = 0U;
  for(auto i = 0U; i < amplicon_count; ++i) {
    if (ampinfo_v[i].graft_cand != no_swarm) {
      ++counter;
    }
  }
  return counter;
}


auto attach_candidates(struct Parameters const & parameters,
                       unsigned int amplicon_count,
                       std::vector<struct ampinfo_s> & ampinfo_v,
                       std::vector<struct swarminfo_s> & swarminfo_v) -> unsigned int
{
  auto const pair_count = count_pairs(amplicon_count, ampinfo_v);

  progress_init("Grafting light swarms on heavy swarms", pair_count);

  /* allocate memory */
  std::vector<struct graft_cand> graft_array(pair_count);

  /* fill in */
  assert(ampinfo_v.size() == amplicon_count);
  auto ticker = 0U;  // refactoring: replace with a transform algorithm
  for(auto i = 0U; i < amplicon_count; ++i) {
    if (ampinfo_v[i].graft_cand == no_swarm) { continue; }
    graft_array[ticker].parent = ampinfo_v[i].graft_cand;
    graft_array[ticker].child = i;  // so two children cannot have the same uint value
    ++ticker;
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
    assert(lhs.child >= rhs.child); // refactoring: child indices are sorted by descending order?
    return lhs.child < rhs.child;
  };

  std::sort(graft_array.begin(), graft_array.end(), compare_grafts);

  /* attach in order */
  auto grafts = 0U;
  auto counter = 1U;
  for(auto const& graft_pair : graft_array) {
      const auto parent = graft_pair.parent;
      const auto child  = graft_pair.child;

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
      progress_update(counter);
      ++counter;
    }
  progress_done(parameters);
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

  for(auto i = 0U; i < variant_count; ++i) {
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

  for(auto i = 0U; i < variant_count; ++i)
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


auto check_heavy_thread(int64_t nth_thread) -> void
{
  static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
  static constexpr auto offset = 4U;
  static constexpr auto nt_per_uint64 = 32U;  // 32 nucleotides can fit in a uint64
  (void) nth_thread;  // refactoring: unused parameter, replace with function overload?

  std::vector<struct var_s> variant_list(multiplier * longestamplicon + offset);
  std::vector<struct var_s> variant_list2(multiplier * (longestamplicon + 1) + offset);

  const std::size_t size =
    sizeof(uint64_t) * ((db_getlongestsequence() + 2 + nt_per_uint64 - 1) / nt_per_uint64);
  std::vector<char> buffer1(size);
  pthread_mutex_lock(&heavy_mutex);
  while ((heavy_amplicon < amplicons) and
         (heavy_progress < heavy_amplicon_count))
    {
      auto const heavy_amplicon_id = heavy_amplicon;
      ++heavy_amplicon;
      assert(heavy_amplicon_id <= std::numeric_limits<std::ptrdiff_t>::max());
      auto const signed_position = static_cast<std::ptrdiff_t>(heavy_amplicon_id);
      auto const & target_amplicon = *std::next(ampinfo, signed_position);
      assert(target_amplicon.swarmid <= std::numeric_limits<std::ptrdiff_t>::max());
      auto const signed_swarmid = static_cast<std::ptrdiff_t>(target_amplicon.swarmid);
      auto const & target_swarm = *std::next(swarminfo, signed_swarmid);
      if (target_swarm.mass >= static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++heavy_progress);  // refactoring: separate operations?
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

  for(auto i = 0U; i < variant_count; ++i) {
    bloomflex_set(bloom, variant_list[i].hash);
  }

  return variant_count;
}


auto mark_light_thread(int64_t nth_thread) -> void
{
  static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
  static constexpr auto offset = 4U;

  (void) nth_thread;  // refactoring: unused?

  std::vector<struct var_s> variant_list(multiplier * longestamplicon + offset);

  pthread_mutex_lock(&light_mutex);
  while (light_progress < light_amplicon_count)
    {
      const auto light_amplicon_id = light_amplicon;
      --light_amplicon;
      assert(light_amplicon_id <= std::numeric_limits<std::ptrdiff_t>::max());
      auto const signed_position = static_cast<std::ptrdiff_t>(light_amplicon_id);
      auto const & target_amplicon = *std::next(ampinfo, signed_position);
      assert(target_amplicon.swarmid <= std::numeric_limits<std::ptrdiff_t>::max());
      auto const signed_swarmid = static_cast<std::ptrdiff_t>(target_amplicon.swarmid);
      auto const & target_swarm = *std::next(swarminfo, signed_swarmid);
      if (target_swarm.mass < static_cast<uint64_t>(opt_boundary))
        {
          progress_update(++light_progress);  // refactoring: separate operations?
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
            if ((opt_no_cluster_breaking) or
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
                    hits_data[hits_count] = amp;
                    ++hits_count;
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

  auto * sequence = db_getsequence(seed);
  const auto seqlen = db_getsequencelen(seed);
  const auto hash = db_gethash(seed);
  const auto variant_count = generate_variants(sequence, seqlen, hash, variant_list);

  // C++17 refactoring:
  // std::for_each_n(variant_list.begin(), variant_count,
  //                 [seed, &hits_data, &hits_count](auto& variant) {
  //                   find_variant_matches(seed, variant, hits_data, hits_count);
  //                 });
  for(auto i = 0U; i < variant_count; ++i) {
    find_variant_matches(seed, variant_list[i], hits_data, hits_count);
  }

  return hits_count;
}


auto network_thread(int64_t nth_thread) -> void
{
  static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
  static constexpr auto offset = 4U;

  (void) nth_thread;  // refactoring: unused?

  std::vector<unsigned int> hits_data(multiplier * longestamplicon + offset + 1);
  std::vector<struct var_s> variant_list(multiplier * longestamplicon + offset + 1);

  pthread_mutex_lock(&network_mutex);
  while (network_amp < amplicons)
    {
      const auto amp = network_amp;
      ++network_amp;
      progress_update(amp);

      pthread_mutex_unlock(&network_mutex);

      const auto hits_count = check_variants(amp, variant_list, hits_data);
      pthread_mutex_lock(&network_mutex);

      assert(amp <= std::numeric_limits<std::ptrdiff_t>::max());
      auto const signed_position = static_cast<std::ptrdiff_t>(amp);
      auto & target_amplicon = *std::next(ampinfo, signed_position);
      target_amplicon.link_start = network_count;
      target_amplicon.link_count = hits_count;

      while (network_count + hits_count > network_v.size()) {
        network_v.reserve(network_v.size() + one_megabyte);
        network_v.resize(network_v.size() + one_megabyte);
      }

      for(auto k = 0U; k < hits_count; ++k) {
        network_v[network_count] = hits_data[k];
        ++network_count;
      }
    }
  pthread_mutex_unlock(&network_mutex);
}


auto process_seed(unsigned int const seed,
                  std::vector<struct ampinfo_s> & ampinfo_v,
                  std::vector<unsigned int> & global_hits_v,
                  unsigned int & global_hits_count) -> void
{
  /* update swarm stats */
  auto const & seed_info = ampinfo_v[seed];

  ++swarmsize;
  swarm_maxgen = std::max(seed_info.generation, swarm_maxgen);
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
        global_hits_alloc += 4UL * one_kilobyte;
      }
      global_hits_v.resize(global_hits_alloc);
      global_hits_data = global_hits_v.data();
    }

  for(auto offset = 0U; offset < link_count; ++offset)
    {
      const auto amp = network_v[link_start + offset];

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


auto compare_amp(const void * void_lhs, const void * void_rhs) -> int
{
  /*
    earlier steps in swarm check that all amplicon sequences are
    unique (strictly dereplicated input data), so distinct amplicons
    with the same sequence are not expected at this stage.
  */
  const auto * lhs = static_cast<const unsigned int*>(void_lhs);
  const auto * rhs = static_cast<const unsigned int*>(void_rhs);
  auto status = 1;  // default is *lhs > *rhs

  assert(*lhs != *rhs);  // '*lhs == *rhs' is not expected at that stage

  // compare amplicon index values (unsigned ints). Amplicon indexes
  // are already sorted by decreasing abundance then by header in
  // db.cc, so smaller indexes should go first. This corresponds to a
  // natural order sorting.
  if (*lhs < *rhs) {
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
  // a network is a cluster with at least two sequences (no singletons)
  progress_init("Dumping network:  ", number_of_networks);

  uint64_t n_processed = 0;  // refactoring: reduce scope (move into the for loop init)
  assert(ampinfo_v.size() == amplicons);
  auto counter = 0ULL;
  for(auto const& amplicon: ampinfo_v) {
    const auto link_start = amplicon.link_start;
    const auto link_count = amplicon.link_count;

    // refactoring: add network tests before replacing with std::sort
    std::qsort(&network_v[link_start],
               link_count,
               sizeof(unsigned int),
               compare_amp);

    // refactoring: std::vector<unsigned int> network_v(network_v.begin() + link_start, network_v.begin() + link_start + link_count);
    for(auto link = 0U; link < link_count; ++link)
      {
        const auto neighbour = network_v[link_start + link];
        fprint_id(parameters.network_file, counter, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        std::fprintf(parameters.network_file, "\t");
        fprint_id(parameters.network_file, neighbour, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        std::fprintf(parameters.network_file, "\n");
        ++n_processed;
      }
    progress_update(n_processed);
    ++counter;
  }
  progress_done(parameters);
}


auto write_swarms_default_format(struct Parameters const & parameters,
                                 std::vector<struct ampinfo_s> & ampinfo_v,
                                 std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  static constexpr char sepchar {' '};
  progress_init("Writing swarms:   ", swarminfo_v.size());

  for(auto i = 0U; i < swarminfo_v.size(); ++i) {
    if (swarminfo_v[i].attached) {
      continue;
    }

    const auto seed = swarminfo_v[i].seed;
    for(auto amp_id = seed; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next) {
      if (amp_id != seed) {
        std::fputc(sepchar, parameters.outfile);
      }
      fprint_id(parameters.outfile, amp_id,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    }
    std::fputc('\n', parameters.outfile);
    progress_update(i + 1);
  }

  progress_done(parameters);
}


auto write_swarms_mothur_format(struct Parameters const & parameters,
                                std::vector<struct ampinfo_s> & ampinfo_v,
                                std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  progress_init("Writing swarms:   ", swarminfo_v.size());

  std::fprintf(parameters.outfile, "swarm_%" PRId64 "\t%" PRIu64,
               parameters.opt_differences, swarmcount_adjusted);

  for(auto i = 0U; i < swarminfo_v.size(); ++i) {
    assert(not swarminfo_v[i].attached);
    if (swarminfo_v[i].attached) {
      continue;
    }

    const auto seed = swarminfo_v[i].seed;
    for(auto amp_id = seed; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next) {
      if (amp_id == seed) {
        std::fputc('\t', parameters.outfile);
      }
      else {
        std::fputc(',', parameters.outfile);
      }
      fprint_id(parameters.outfile, amp_id,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    }
    progress_update(i + 1);
  }

  std::fputc('\n', parameters.outfile);

  progress_done(parameters);
}


auto write_swarms_uclust_format(struct Parameters const & parameters,
                                std::vector<struct ampinfo_s> & ampinfo_v,
                                std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  static constexpr double one_hundred = 100.0;
  auto cluster_no = 0U;
  const auto score_matrix_63 = create_score_matrix<int64_t>(parameters.penalty_mismatch);
  std::vector<unsigned char> directions(1UL * longestamplicon * longestamplicon);
  std::vector<uint64_t> hearray(2UL * longestamplicon);
  std::vector<char> raw_alignment;
  std::string cigar_string;
  raw_alignment.reserve(2UL * longestamplicon);
  cigar_string.reserve(2UL * longestamplicon);

  progress_init("Writing UCLUST:   ", swarminfo_v.size());

  auto counter = 0U;
  for (auto const & swarm_info : swarminfo_v) {
    if (swarm_info.attached) {
      continue;
    }

    const auto seed = swarm_info.seed;

    auto const & seed_info = ampinfo_v[seed];

    std::fprintf(parameters.uclustfile, "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                 cluster_no,
                 swarm_info.size);
    fprint_id(parameters.uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    std::fprintf(parameters.uclustfile, "\t*\n");

    std::fprintf(parameters.uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                 cluster_no,
                 db_getsequencelen(seed));
    fprint_id(parameters.uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    std::fprintf(parameters.uclustfile, "\t*\n");

    for(auto amp_id = seed_info.next; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next)
      {
        auto * dseq = db_getsequence(amp_id);
        const auto dlen = db_getsequencelen(amp_id);  // refactoring: as a struct Sequence{ptr, length}
        auto * qseq = db_getsequence(seed);  // refactoring: can be moved outside of this loop!
        const auto qlen = db_getsequencelen(seed);

        uint64_t nwdiff = 0;  // refactoring: nw() -> uint64_t?

        nw(dseq, dlen, qseq, qlen,
           score_matrix_63, static_cast<unsigned long int>(parameters.penalty_gapopen),
           static_cast<unsigned long int>(parameters.penalty_gapextend),
           nwdiff, directions, hearray, raw_alignment);

        // backtracking produces a reversed alignment (starting from the end)
        std::reverse(raw_alignment.begin(), raw_alignment.end());
        compress_alignment_to_cigar(raw_alignment, cigar_string);

        // loosing precision when converting raw_alignment.size() and
        // nwdiff to double is not an issue, no need to add assertions
        const auto nwalignmentlength = static_cast<double>(raw_alignment.size());
        const auto differences = static_cast<double>(nwdiff);
        const double percentid = one_hundred * (nwalignmentlength - differences) / nwalignmentlength;

        std::fprintf(parameters.uclustfile,
                     "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                     cluster_no,
                     db_getsequencelen(amp_id),
                     percentid,
                     nwdiff > 0 ? cigar_string.data() : "=");

        fprint_id(parameters.uclustfile, amp_id, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        std::fprintf(parameters.uclustfile, "\t");
        fprint_id(parameters.uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        std::fprintf(parameters.uclustfile, "\n");

        raw_alignment.clear();
        cigar_string.clear();
      }

    ++cluster_no;
    progress_update(counter);
    ++counter;
  }
  progress_done(parameters);
}


auto write_representative_sequences(struct Parameters const & parameters,
                                    std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  progress_init("Writing seeds:    ", swarminfo_v.size());

  std::vector<unsigned int> sorter(swarminfo_v.size());
  std::iota(sorter.begin(), sorter.end(), 0);

  auto compare_mass_and_headers = [&swarminfo_v](unsigned int const lhs,
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
    // assert(status != 0); // all headers are unique
    return status < 0;
  };

  std::sort(sorter.begin(), sorter.end(), compare_mass_and_headers);

  auto counter = 1U;
  for (const auto index : sorter) {
    const auto & a_swarm = swarminfo_v[index];
    if (a_swarm.attached) {
      continue;
    }
    const auto seed = a_swarm.seed;
    const auto mass = a_swarm.mass;
    std::fprintf(parameters.seeds_file, ">");
    fprint_id_with_new_abundance(parameters.seeds_file, seed, mass,
                                 parameters.opt_usearch_abundance);
    std::fprintf(parameters.seeds_file, "\n");
    db_fprintseq(parameters.seeds_file, seed);
    progress_update(counter);
    ++counter;
  }

  progress_done(parameters);
}


auto write_structure_file(struct Parameters const & parameters,
                          std::vector<struct ampinfo_s> & ampinfo_v,
                          std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  auto cluster_no = 0U;

  progress_init("Writing structure:", swarminfo_v.size());

  for(auto swarmid = 0U; swarmid < swarminfo_v.size(); ++swarmid)
    {
      if (swarminfo_v[swarmid].attached) {
        continue;
      }
      const auto seed = swarminfo_v[swarmid].seed;

      auto const & seed_info = ampinfo_v[seed];

      for(auto amp_id = seed_info.next; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next)
        {
          const auto graft_parent = ampinfo_v[amp_id].graft_cand;
          if (graft_parent != no_swarm)
            {
              fprint_id_noabundance(parameters.internal_structure_file,
                                    graft_parent, parameters.opt_usearch_abundance);
              std::fprintf(parameters.internal_structure_file, "\t");
              fprint_id_noabundance(parameters.internal_structure_file, amp_id, parameters.opt_usearch_abundance);
              std::fprintf(parameters.internal_structure_file,
                           "\t%d\t%u\t%u\n",
                           2,
                           cluster_no + 1,
                           ampinfo_v[graft_parent].generation + 1);
            }

          const auto parent = ampinfo_v[amp_id].parent;
          if (parent != no_swarm)
            {
              fprint_id_noabundance(parameters.internal_structure_file, parent, parameters.opt_usearch_abundance);
              std::fprintf(parameters.internal_structure_file, "\t");
              fprint_id_noabundance(parameters.internal_structure_file, amp_id, parameters.opt_usearch_abundance);
              std::fprintf(parameters.internal_structure_file,
                           "\t%u\t%u\t%u\n",
                           1U,
                           cluster_no + 1,
                           ampinfo_v[amp_id].generation);
            }
        }

      ++cluster_no;
      progress_update(swarmid);
    }
  progress_done(parameters);
}


auto write_stats_file(struct Parameters const & parameters,
                      std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  progress_init("Writing stats:    ", swarminfo_v.size());

  auto counter = 0U;
  for (auto const & swarm_info : swarminfo_v) {
    assert(not swarm_info.attached);
    if (swarm_info.attached) {
      continue;
    }
    std::fprintf(parameters.statsfile, "%u\t%" PRIu64 "\t", swarm_info.size, swarm_info.mass);
    fprint_id_noabundance(parameters.statsfile, swarm_info.seed, parameters.opt_usearch_abundance);
    std::fprintf(parameters.statsfile, "\t%" PRIu64 "\t%u\t%u\t%u\n",
                 db_getabundance(swarm_info.seed),
                 swarm_info.singletons, swarm_info.maxgen, swarm_info.maxgen);
    progress_update(counter);
    ++counter;
  }
  progress_done(parameters);
}


auto output_results(struct Parameters const & parameters,
                    std::vector<struct ampinfo_s> & ampinfo_v,
                    std::vector<struct swarminfo_s> & swarminfo_v) -> void {
  /* dump swarms */
  if (parameters.opt_mothur) {
    write_swarms_mothur_format(parameters, ampinfo_v, swarminfo_v);
  }
  else {
    write_swarms_default_format(parameters, ampinfo_v, swarminfo_v);
  }

  /* dump seeds in fasta format with sum of abundances */
  if (not parameters.opt_seeds.empty()) {
    write_representative_sequences(parameters, swarminfo_v);
  }

  /* output internal structure */
  if (not parameters.opt_internal_structure.empty()) {
    write_structure_file(parameters, ampinfo_v, swarminfo_v);
  }

  /* output swarms in uclust format */
  if (not parameters.opt_uclust_file.empty()) {
    write_swarms_uclust_format(parameters, ampinfo_v, swarminfo_v);
  }

  /* output statistics to file */
  if (not parameters.opt_statistics_file.empty()) {
    write_stats_file(parameters, swarminfo_v);
  }
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
  static constexpr auto multiplier = 7U;
  static constexpr auto offset = 4U;
  const auto global_hits_alloc = multiplier * longestamplicon + offset + 1;
  std::vector<unsigned int> global_hits_v(global_hits_alloc);
  global_hits_data = global_hits_v.data();

  /* compute hash for all amplicons and store them in a hash table */
  std::vector<unsigned char> hash_occupied_v;
  std::vector<uint64_t> hash_values_v;
  std::vector<unsigned int> hash_data_v;
  const auto hashtablesize = hash_alloc(amplicons,
                                        hash_occupied_v,
                                        hash_values_v,
                                        hash_data_v);
  struct bloom_s bloom_filter;
  bloom_a = bloom_init(hashtablesize, bloom_filter);

  duplicates_found = 0;

  progress_init("Hashing sequences:", amplicons);
  for(auto k = 0U; k < amplicons; ++k)
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

  progress_done(parameters);

  /* for all amplicons, generate list of matching amplicons */
  network_v.resize(one_megabyte);

  network_count = 0;

  pthread_mutex_init(&network_mutex, nullptr);
  network_amp = 0;
  progress_init("Building network: ", amplicons);
  {
    assert(parameters.opt_threads <= std::numeric_limits<int>::max());
    // refactoring C++14: use std::make_unique
    std::unique_ptr<ThreadRunner> network_tr (new ThreadRunner(static_cast<int>(parameters.opt_threads), network_thread));
    network_tr->run();
  }
  pthread_mutex_destroy(&network_mutex);

  progress_done(parameters);


  /* dump network to file */
  if (not parameters.opt_network_file.empty()) {
    write_network_file(network_count, parameters, ampinfo_v);
  }


  /* for each non-swarmed amplicon look for subseeds ... */

  auto swarmcount = 0U;  // refactoring: find a way to know swarmcount in advance?
  progress_init("Clustering:       ", amplicons);

  for(auto seed = 0U; seed < amplicons; ++seed)
    {
      auto & seed_info = ampinfo_v[seed];

      if (seed_info.swarmid == no_swarm)
        {
          /* start a new swarm with a new initial seed */

          seed_info.swarmid = swarmcount;
          seed_info.generation = 0;
          seed_info.parent = no_swarm;
          seed_info.next = no_swarm;

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
          for(auto i = 0U; i < global_hits_count; ++i) {
            add_amp_to_swarm(global_hits_v[i], ampinfo_v);
          }

          /* find later generation matches */
          auto subseed = seed_info.next;
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
              for(auto i = 0U; i < global_hits_count; ++i) {
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

          if (swarmcount >= swarminfo_v.size())
            {
              /* allocate memory for more swarms... */
              // note: capacity doubles, as usual
              // 1,024 times struct size (so at least 40,960 new bytes reserved)
              swarminfo_v.resize(swarminfo_v.size() + one_kilobyte);
              swarminfo = swarminfo_v.data();
            }

          auto & swarm_info = swarminfo_v[swarmcount];

          swarm_info.seed = seed;
          swarm_info.size = swarmsize;
          swarm_info.mass = abundance_sum;
          swarm_info.sumlen = swarm_sumlen;
          swarm_info.singletons = singletons;
          swarm_info.maxgen = swarm_maxgen;
          swarm_info.last = current_swarm_tail;
          swarm_info.attached = false;

          /* update overall stats */
          largest = std::max(swarmsize, largest);
          maxgen = std::max(swarm_maxgen, maxgen);

          ++swarmcount;
        }
      progress_update(seed + 1);
    }
  progress_done(parameters);

  global_hits_data = nullptr;

  network_v.clear();
  network_v.shrink_to_fit();

  swarmcount_adjusted = swarmcount;

  /* fastidious */

  if (parameters.opt_fastidious)
    {
      std::fprintf(parameters.logfile, "\n");
      std::fprintf(parameters.logfile, "Results before fastidious processing:\n");
      std::fprintf(parameters.logfile, "Number of swarms:  %u\n", swarmcount);
      std::fprintf(parameters.logfile, "Largest swarm:     %u\n", largest);
      std::fprintf(parameters.logfile, "\n");

      uint64_t small_clusters = 0;
      uint64_t amplicons_in_small_clusters = 0;
      uint64_t nucleotides_in_small_clusters = 0;

      // refactoring: move to function that returns a struct cluster_stats
      progress_init("Counting amplicons in heavy and light swarms",
                    swarmcount);

      for(auto i = 0ULL; i < swarmcount; ++i)
        {
          auto const & swarm_info = swarminfo_v[i];
          if (swarm_info.mass < static_cast<uint64_t>(opt_boundary))
            {
              amplicons_in_small_clusters += swarm_info.size;
              nucleotides_in_small_clusters += swarm_info.sumlen;
              ++small_clusters;
            }
          progress_update(i + 1);
        }
      progress_done(parameters);

      const uint64_t amplicons_in_large_clusters = amplicons - amplicons_in_small_clusters;
      const uint64_t large_clusters = swarmcount - small_clusters;

      std::fprintf(parameters.logfile, "Heavy swarms: %" PRIu64 ", with %" PRIu64 " amplicons\n",
                   large_clusters, amplicons_in_large_clusters);
      std::fprintf(parameters.logfile, "Light swarms: %" PRIu64 ", with %" PRIu64 " amplicons\n",
                   small_clusters, amplicons_in_small_clusters);
      std::fprintf(parameters.logfile, "Total length of amplicons in light swarms: %" PRIu64 "\n",
                   nucleotides_in_small_clusters);

      if ((small_clusters == 0) or (large_clusters == 0))
        {
          std::fprintf(parameters.logfile, "Only light or heavy swarms found - "
                       "no need for further analysis.\n");
        }
      else
        {
          /* m: total size of Bloom filter in bits */
          /* k: number of hash functions (n_hash_functions) */
          /* n: number of entries in the bloom filter */
          /* here: k=11 and m/n=18, that is 16 bits/entry */

          static constexpr auto microvariants = 7U;
          static constexpr auto n_bits_in_a_byte = 8U;
          static constexpr double hash_functions_per_bit {4.0 / 10};
          static constexpr double natural_log_of_2 {0.693147181};  // C++26 refactoring: std::log(2.0)
          static_assert(hash_functions_per_bit <= natural_log_of_2, "upper limit is log(2)");
          assert(parameters.opt_bloom_bits <= std::numeric_limits<unsigned int>::max());
          assert(parameters.opt_bloom_bits <= 64);  // larger than expected
          assert(parameters.opt_bloom_bits >= 2);  // smaller than expected
          auto bits = static_cast<uint64_t>(parameters.opt_bloom_bits);
          auto bits_uint = static_cast<unsigned int>(parameters.opt_bloom_bits);  // avoid risky conversion warning: uint64 to double

          // int64_t n_hash_functions = int(bits * std::log(2.0));    /* 16 bits -> 11 hash functions */
          // auto n_hash_functions = unsigned int(hash_functions_per_bit * bits); /* 6 */
          auto n_hash_functions = std::max(static_cast<unsigned int>(hash_functions_per_bit * bits_uint), 1U);

          uint64_t bloom_length_in_bits = nucleotides_in_small_clusters * microvariants * bits;

          const uint64_t memtotal = arch_get_memtotal();
          const uint64_t memused = arch_get_memused();

          if (parameters.opt_ceiling != 0)
            {
              if (static_cast<uint64_t>(parameters.opt_ceiling) * one_megabyte < memused)
                {
                  fatal(error_prefix, "Memory ceiling for Bloom filter is too low.");
                }
              assert(memused < one_megabyte * static_cast<uint64_t>(parameters.opt_ceiling));
              const uint64_t memrest
                = one_megabyte * static_cast<uint64_t>(parameters.opt_ceiling) - memused;
              auto const new_bits = n_bits_in_a_byte * memrest / (microvariants * nucleotides_in_small_clusters);
              if (new_bits < bits)
                {
                  if (new_bits < 2) {
                    fatal(error_prefix, "Insufficient memory remaining for Bloom filter.");
                  }
                  std::fprintf(parameters.logfile, "Reducing memory used for Bloom filter due to --ceiling option.\n");
                  bits = new_bits;
                  bits_uint = static_cast<unsigned int>(new_bits);
                  n_hash_functions = std::max(static_cast<unsigned int>(hash_functions_per_bit * bits_uint), 1U);
                  bloom_length_in_bits = nucleotides_in_small_clusters * microvariants * bits;
                }
            }

          static constexpr uint64_t min_bloom_length_in_bits {64};  // at least 64 bits
          bloom_length_in_bits = std::max(bloom_length_in_bits, min_bloom_length_in_bits);

          if (memused + bloom_length_in_bits / n_bits_in_a_byte > memtotal)
            {
              std::fprintf(parameters.logfile, "WARNING: Memory usage will probably exceed total amount of memory available.\n");
              std::fprintf(parameters.logfile, "Try to reduce memory footprint using the --bloom-bits or --ceiling options.\n");
            }

          std::fprintf(parameters.logfile,
                       "Bloom filter: bits=%" PRIu64 ", m=%" PRIu64 ", k=%u, size=%.1fMB\n",
                       bits, bloom_length_in_bits, n_hash_functions, static_cast<double>(bloom_length_in_bits) / (n_bits_in_a_byte * one_megabyte));


          // bloom_length is in bits (divide by 8 to get bytes)
          // bloom_length is guaranteed to be at least 64 (see code above)
          assert(bloom_length_in_bits != 0);  // safeguard for future changes
          assert(bloom_length_in_bits >= 64);
          const uint64_t n_bytes = ((bloom_length_in_bits - 1) / n_bits_in_a_byte) + 1;
          struct bloomflex_s bloomflex_filter;
          bloom_f = bloomflex_init(n_bytes, n_hash_functions, bloomflex_filter);


          /* Empty the old hash and bloom filter
             before we reinsert only the light swarm amplicons */

          std::fill(hash_occupied_v.begin(), hash_occupied_v.end(), 0U);
          bloom_zap(bloom_filter);

          progress_init("Adding light swarm amplicons to Bloom filter",
                        amplicons_in_small_clusters);

          /* process amplicons in order from least to most abundant */
          /* but stop when all amplicons in small clusters are processed */

          light_variants = 0;

          pthread_mutex_init(&light_mutex, nullptr);
          light_progress = 0;
          light_amplicon_count = amplicons_in_small_clusters;
          light_amplicon = amplicons - 1;
          {
            assert(parameters.opt_threads <= std::numeric_limits<int>::max());
            // refactoring C++14: use std::make_unique
            std::unique_ptr<ThreadRunner> light_tr (new ThreadRunner(static_cast<int>(parameters.opt_threads), mark_light_thread));
            light_tr->run();
          }
          pthread_mutex_destroy(&light_mutex);

          progress_done(parameters);

          std::fprintf(parameters.logfile,
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
            assert(parameters.opt_threads <= std::numeric_limits<int>::max());
            // refactoring C++14: use std::make_unique
            std::unique_ptr<ThreadRunner> heavy_tr (new ThreadRunner(static_cast<int>(parameters.opt_threads), check_heavy_thread));
            heavy_tr->run();
          }

          pthread_mutex_destroy(&heavy_mutex);

          progress_done(parameters);

          bloomflex_exit(bloomflex_filter);

          pthread_mutex_destroy(&graft_mutex);

          std::fprintf(parameters.logfile, "Heavy variants: %" PRIu64 "\n", heavy_variants);
          std::fprintf(parameters.logfile, "Got %" PRId64 " graft candidates\n", graft_candidates);
          const unsigned int grafts = attach_candidates(parameters, amplicons, ampinfo_v, swarminfo_v);
          std::fprintf(parameters.logfile, "Made %u grafts\n", grafts);
          std::fprintf(parameters.logfile, "\n");
        }
    }

  // refactoring: trim vectors (remove allocated unused elements)
  // could it be done before the fastidious phase?
  swarminfo_v.resize(swarmcount);  // swarminfo_v's capacity can be twice too much
  swarminfo_v.shrink_to_fit();

  output_results(parameters, ampinfo_v, swarminfo_v);

  std::fprintf(parameters.logfile, "\n");
  std::fprintf(parameters.logfile, "Number of swarms:  %" PRIu64 "\n", swarmcount_adjusted);
  std::fprintf(parameters.logfile, "Largest swarm:     %u\n", largest);
  std::fprintf(parameters.logfile, "Max generations:   %u\n", maxgen);

  bloom_exit(bloom_a);
  hash_free();

  swarminfo = nullptr;
  ampinfo = nullptr;
}
