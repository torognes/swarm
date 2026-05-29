/*
    SWARM

    Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

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
#include "utils/system_memory.h"
#include "utils/bloom.h"
#include "db.h"
#include "utils/hashtable.h"
#include "utils/nw_aligner.h"
#include "variants.h"
#include "utils/fatal.h"
#include "utils/make_unique.h"
#include "utils/nt_codec.h"
#include "utils/progress.h"
#include "utils/threads.h"
#include "utils/view.h"
#include <algorithm>  // std::sort(), std::max()
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc(), size_t
#include <iterator>  // std::next
#include <limits>  // unsigned int max
#include <memory>  // unique pointer
#include <mutex>  // std::lock_guard, std::unique_lock
#include <numeric>  // std::iota
#include <string>
#include <vector>


namespace {

  constexpr unsigned int one_kilobyte {1U << 10U};  // 1,024 bytes
  constexpr unsigned int one_megabyte {one_kilobyte * one_kilobyte};
  constexpr unsigned int no_swarm {std::numeric_limits<unsigned int>::max()};


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

  struct graft_cand
  {
    unsigned int parent;
    unsigned int child;
  };

  /* Information about potential grafts */
  struct Graft_state
  {
    std::mutex mutex;
    int64_t candidates {0};
  };

  /* overall statistics, accumulated across all swarms */
  struct Overall_stats
  {
    uint64_t swarmcount_adjusted {0};
    unsigned int maxgen {0};
    unsigned int largest {0};
  };

  /* per-swarm statistics, reset before each new swarm and copied into
     the corresponding swarminfo_s entry once the swarm is closed */
  struct Active_swarm_stats
  {
    uint64_t abundance_sum {0};  /* = mass */
    uint64_t sumlen {0};
    unsigned int tail {0};
    unsigned int size {0};
    unsigned int maxgen {0};
    unsigned int singletons {0};
  };

  Overall_stats overall_stats {};
  Active_swarm_stats current_swarm {};

  unsigned int * global_hits_data {nullptr};

  unsigned int amplicons {0};

  struct Heavy_state
  {
    std::mutex mutex;
    uint64_t variants {0};
    uint64_t progress {0};
    uint64_t amplicon_count {0};
    unsigned int amplicon {0};
  };

  struct Light_state
  {
    std::mutex mutex;
    uint64_t variants {0};
    uint64_t progress {0};
    uint64_t amplicon_count {0};
    unsigned int amplicon {0};
  };

  struct Network_state
  {
    std::mutex mutex;
    unsigned int amp {0};
    unsigned int count {0};
    std::vector<unsigned int> network_v;
  };

  struct Cluster_stats
  {
    uint64_t small_clusters {0};
    uint64_t large_clusters {0};
    uint64_t amplicons_in_small_clusters {0};
    uint64_t amplicons_in_large_clusters {0};
    uint64_t nucleotides_in_small_clusters {0};
  };

  struct Bloom_geometry
  {
    uint64_t n_bytes {0};
    unsigned int n_hash_functions {0};
  };

  /* Bloom filter shape used for the per-amplicon hashtable + bloom_a
     in both the d=1 phase and the fastidious phase. */
  constexpr unsigned int amplicon_pattern_shift {10};
  constexpr unsigned int amplicon_n_hash_functions {8};


  inline auto hash_insert(Data const & data,
                          Hashtable & hash_table,
                          BloomFilter & bloom_a,
                          unsigned int const amp) -> void {
    /* find the first empty bucket */
    const auto hash = data.sequence_hash(amp);
    auto index = hash_table.getindex(hash);
    while (hash_table.is_occupied(index)) {
      index = hash_table.getnextindex(index);
    }

    hash_table.set_occupied(index);
    hash_table.set_value(index, hash);
    hash_table.set_data(index, amp);
    bloom_a.set(hash);
  }


  /******************** FASTIDIOUS START ********************/


  auto attach(unsigned int seed, unsigned int amp,
              std::vector<struct ampinfo_s> & ampinfo_v,
              std::vector<struct swarminfo_s> & swarminfo_v) -> void
  {
    /* graft light swarm (amp) on heavy swarm (seed) */

    auto & heavy_swarm = swarminfo_v[ampinfo_v[seed].swarmid];
    auto & light_swarm = swarminfo_v[ampinfo_v[amp].swarmid];

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
    overall_stats.largest = std::max(heavy_swarm.size, overall_stats.largest);

    --overall_stats.swarmcount_adjusted;
  }


  auto add_graft_candidate(std::vector<struct ampinfo_s> & ampinfo_v,
                           unsigned int seed, unsigned int amp,
                           struct Graft_state & graft_state) -> void
  {
    std::lock_guard<std::mutex> const lock(graft_state.mutex);
    ++graft_state.candidates;
    auto & amplicon = ampinfo_v[amp];
    // if there is no heavy candidate to graft amp, or if seed is
    // earlier in the sorting order, then we change the attachment to
    // seed
    if ((amplicon.graft_cand == no_swarm) or (amplicon.graft_cand > seed)) {
      amplicon.graft_cand = seed;
    }
  }


  // C++17 refactoring: replace with std::count_if()
  auto count_pairs(std::vector<struct ampinfo_s> const & ampinfo_v) -> unsigned int {
    auto counter = 0U;
    for (auto const & info : ampinfo_v) {
      if (info.graft_cand != no_swarm) {
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
    auto const pair_count = count_pairs(ampinfo_v);

    Progress progress("Grafting light swarms on heavy swarms", pair_count, parameters);

    /* allocate memory */
    std::vector<struct graft_cand> graft_array(pair_count);

    /* fill in */
    assert(ampinfo_v.size() == amplicon_count);
    auto ticker = 0U;  // refactoring: replace with a transform algorithm
    for (auto i = 0U; i < amplicon_count; ++i) {
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
    for (auto const& graft_pair : graft_array) {
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
      progress.update();
    }
    progress.done();
    return grafts;
  }


  auto hash_check_attach(Data const & data,
                         std::vector<struct ampinfo_s> & ampinfo_v,
                         Hashtable const & hash_table,
                         Sequence const & seed_seq,
                         struct var_s const & var,
                         unsigned int seed,
                         struct Graft_state & graft_state) -> bool
  {
    /* seed is the original large swarm seed */

    /* compute hash and corresponding hash table index */
    const auto hash = var.hash;
    auto index = hash_table.getindex(hash);

    /* find matching buckets */

    while (hash_table.is_occupied(index))
      {
        if (hash_table.compare_value(index, hash))
          {
            /* check that mass is below threshold */
            const auto amp = hash_table.get_data(index);

            /* make absolutely sure sequences are identical */
            auto const amp_seq = data.sequence_view(amp);
            if (check_variant(seed_seq, var, amp_seq))
              {
                add_graft_candidate(ampinfo_v, seed, amp, graft_state);
                return true;
              }
          }
        index = hash_table.getnextindex(index);
      }
    return false;
  }


  inline auto check_heavy_var_2(Data const & data,
                                std::vector<struct ampinfo_s> & ampinfo_v,
                                Hashtable const & hash_table,
                                BloomFilter const & bloom_a,
                                Sequence const & seq,
                                unsigned int seed,
                                std::vector<struct var_s>& variant_list,
                                struct Graft_state & graft_state) -> uint64_t
  {
    /* Check second generation microvariants of the heavy swarm amplicons
       and see if any of them are identical to a light swarm amplicon. */

    uint64_t matches = 0;

    const auto hash = data.zobrist().hash(seq);
    const auto variant_count = generate_variants(data.zobrist(), seq, hash, variant_list);

    // variant_list is pre-sized to an upper bound; only the first
    // variant_count entries are valid for this call.
    for (auto const & var : View<var_s>{variant_list.data(), variant_count}) {
      if (bloom_a.get(var.hash) and
          hash_check_attach(data, ampinfo_v, hash_table, seq, var, seed, graft_state)) {
        ++matches;
      }
    }

    return matches;
  }


  auto check_heavy_var(Data const & data,
                       std::vector<struct ampinfo_s> & ampinfo_v,
                       Hashtable const & hash_table,
                       BloomFilter const & bloom_a,
                       BloomFilter const & bloom_f,
                       std::vector<char>& varseq,
                       unsigned int seed,
                       uint64_t & number_of_matches,
                       uint64_t & number_of_variants,
                       std::vector<struct var_s>& variant_list,
                       std::vector<struct var_s>& variant_list2,
                       struct Graft_state & graft_state) -> void
  {
    /*
      bloom_f is the fastidious bloom filter in which to check the variants
      bloom_a is the amplicon bloom filter, used in check_heavy_var_2
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

    auto const seed_seq = data.sequence_view(seed);
    const auto hash = data.sequence_hash(seed);
    const auto variant_count = generate_variants(data.zobrist(), seed_seq, hash, variant_list);

    for (auto i = 0U; i < variant_count; ++i)
      {
        struct var_s const & var = variant_list[i];
        if (bloom_f.get(var.hash))
          {
            auto varlen = 0U;
            generate_variant_sequence(seed_seq, var, varseq, varlen);
            auto const var_seq = Sequence{View<char>{varseq.data(), nt_bytelength(varlen)}, varlen};
            matches += check_heavy_var_2(data, ampinfo_v, hash_table,
                                         bloom_a,
                                         var_seq,
                                         seed,
                                         variant_list2,
                                         graft_state);
          }
      }

    number_of_matches = matches;
    number_of_variants = variant_count;
  }


  auto check_heavy_thread(struct Parameters const & parameters,
                          Data const & data,
                          std::vector<struct ampinfo_s> & ampinfo_v,
                          std::vector<struct swarminfo_s> const & swarminfo_v,
                          Hashtable const & hash_table,
                          BloomFilter const & bloom_a,
                          BloomFilter const & bloom_f,
                          uint64_t nth_thread,
                          struct Heavy_state & heavy_state,
                          struct Graft_state & graft_state,
                          Progress & progress) -> void
  {
    static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
    static constexpr auto offset = 4U;
    static constexpr auto nt_per_uint64 = 32U;  // 32 nucleotides can fit in a uint64
    (void) nth_thread;  // refactoring: unused parameter, replace with function overload?

    std::vector<struct var_s> variant_list((multiplier * data.longest_sequence()) + offset);
    std::vector<struct var_s> variant_list2((multiplier * (data.longest_sequence() + 1)) + offset);

    const std::size_t size =
      sizeof(uint64_t) * ((data.longest_sequence() + 2 + nt_per_uint64 - 1) / nt_per_uint64);
    std::vector<char> buffer1(size);
    std::unique_lock<std::mutex> lock(heavy_state.mutex);
    while ((heavy_state.amplicon < amplicons) and
           (heavy_state.progress < heavy_state.amplicon_count))
      {
        auto const heavy_amplicon_id = heavy_state.amplicon;
        ++heavy_state.amplicon;
        auto const & target_amplicon = ampinfo_v[heavy_amplicon_id];
        auto const & target_swarm = swarminfo_v[target_amplicon.swarmid];
        if (target_swarm.mass >= static_cast<uint64_t>(parameters.opt_boundary))
          {
            progress.update(++heavy_state.progress);  // refactoring: separate operations?
            lock.unlock();
            uint64_t number_of_matches {0};
            uint64_t number_of_variants {0};
            check_heavy_var(data, ampinfo_v, hash_table, bloom_a, bloom_f, buffer1, heavy_amplicon_id,
                            number_of_matches, number_of_variants,
                            variant_list, variant_list2,
                            graft_state);
            lock.lock();
            heavy_state.variants += number_of_variants;
          }
      }
  }


  auto mark_light_var(Data const & data,
                      Hashtable & hash_table,
                      BloomFilter & bloom_a,
                      BloomFilter & bloom_f,
                      unsigned int seed,
                      std::vector<struct var_s>& variant_list) -> uint64_t
  {
    /*
      add all microvariants of seed to Bloom filter

      bloom_f is the fastidious BloomFilter in which to enter the variants
      bloom_a is the amplicon BloomFilter, updated by hash_insert
      seed is the original seed
    */

    hash_insert(data, hash_table, bloom_a, seed);

    auto const seed_seq = data.sequence_view(seed);
    const auto hash = data.sequence_hash(seed);
    const auto variant_count = generate_variants(data.zobrist(), seed_seq, hash, variant_list);

    for (auto i = 0U; i < variant_count; ++i) {
      bloom_f.set(variant_list[i].hash);
    }

    return variant_count;
  }


  auto mark_light_thread(struct Parameters const & parameters,
                         Data const & data,
                         std::vector<struct ampinfo_s> const & ampinfo_v,
                         std::vector<struct swarminfo_s> const & swarminfo_v,
                         Hashtable & hash_table,
                         BloomFilter & bloom_a,
                         BloomFilter & bloom_f,
                         uint64_t nth_thread,
                         struct Light_state & state,
                         Progress & progress) -> void
  {
    static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
    static constexpr auto offset = 4U;

    (void) nth_thread;  // refactoring: unused?

    std::vector<struct var_s> variant_list((multiplier * data.longest_sequence()) + offset);

    std::unique_lock<std::mutex> lock(state.mutex);
    while (state.progress < state.amplicon_count)
      {
        const auto light_amplicon_id = state.amplicon;
        --state.amplicon;
        auto const & target_amplicon = ampinfo_v[light_amplicon_id];
        auto const & target_swarm = swarminfo_v[target_amplicon.swarmid];
        if (target_swarm.mass < static_cast<uint64_t>(parameters.opt_boundary))
          {
            progress.update(++state.progress);  // refactoring: separate operations?
            lock.unlock();
            const auto variant_count = mark_light_var(data, hash_table, bloom_a, bloom_f,
                                                      light_amplicon_id,
                                                      variant_list);
            lock.lock();
            state.variants += variant_count;
          }
      }
  }


  /******************** FASTIDIOUS END ********************/


  inline auto find_variant_matches(struct Parameters const & parameters,
                                   Data const & data,
                                   Hashtable const & hash_table,
                                   BloomFilter const & bloom_a,
                                   unsigned int seed,
                                   struct var_s const & var,
                                   std::vector<unsigned int>& hits_data,
                                   unsigned int & hits_count) -> void
  {
    if (not bloom_a.get(var.hash)) {
      return;
    }

    /* compute hash and corresponding hash table index */

    auto index = hash_table.getindex(var.hash);

    /* find matching buckets */

    while (hash_table.is_occupied(index))
      {
        if (hash_table.compare_value(index, var.hash))
          {
            const auto amp = hash_table.get_data(index);

            /* avoid self */
            if (seed != amp) {
              if ((parameters.opt_no_cluster_breaking) or
                  (data.abundance(seed) >= data.abundance(amp)))
                {
                  auto const seed_seq = data.sequence_view(seed);
                  auto const amp_seq = data.sequence_view(amp);

                  if (check_variant(seed_seq, var, amp_seq))
                    {
                      hits_data[hits_count] = amp;
                      ++hits_count;
                      break;
                    }
                }
            }
          }
        index = hash_table.getnextindex(index);
      }
  }


  auto check_variants(struct Parameters const & parameters,
                      Data const & data,
                      Hashtable const & hash_table,
                      BloomFilter const & bloom_a,
                      unsigned int seed,
                      std::vector<struct var_s> & variant_list,
                      std::vector<unsigned int>& hits_data) -> unsigned int
  {
    auto hits_count = 0U;

    auto const seed_seq = data.sequence_view(seed);
    const auto hash = data.sequence_hash(seed);
    const auto variant_count = generate_variants(data.zobrist(), seed_seq, hash, variant_list);

    // C++17 refactoring:
    // std::for_each_n(variant_list.begin(), variant_count,
    //                 [seed, &hits_data, &hits_count](auto& variant) {
    //                   find_variant_matches(parameters, hash_table, bloom_a, seed, variant, hits_data, hits_count);
    //                 });
    for (auto i = 0U; i < variant_count; ++i) {
      find_variant_matches(parameters, data, hash_table, bloom_a, seed, variant_list[i], hits_data, hits_count);
    }

    return hits_count;
  }


  auto network_thread(struct Parameters const & parameters,
                      Data const & data,
                      std::vector<struct ampinfo_s> & ampinfo_v,
                      Hashtable const & hash_table,
                      BloomFilter const & bloom_a,
                      uint64_t nth_thread,
                      struct Network_state & state,
                      Progress & progress) -> void
  {
    static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
    static constexpr auto offset = 4U;
    std::size_t const n_items = (multiplier * data.longest_sequence()) + offset + 1;

    (void) nth_thread;  // refactoring: unused?

    std::vector<unsigned int> hits_data(n_items);
    std::vector<struct var_s> variant_list(n_items);

    std::unique_lock<std::mutex> lock(state.mutex);
    while (state.amp < amplicons)
      {
        const auto amp = state.amp;
        ++state.amp;
        progress.update(amp);

        lock.unlock();

        const auto hits_count = check_variants(parameters, data, hash_table, bloom_a, amp, variant_list, hits_data);
        lock.lock();

        auto & target_amplicon = ampinfo_v[amp];
        target_amplicon.link_start = state.count;
        target_amplicon.link_count = hits_count;

        while (state.count + hits_count > state.network_v.size()) {
          state.network_v.reserve(state.network_v.size() + one_megabyte);
          state.network_v.resize(state.network_v.size() + one_megabyte);
        }

        for (auto k = 0U; k < hits_count; ++k) {
          state.network_v[state.count] = hits_data[k];
          ++state.count;
        }
      }
  }


  auto process_seed(Data const & data,
                    unsigned int const seed,
                    std::vector<struct ampinfo_s> & ampinfo_v,
                    std::vector<unsigned int> const & network_v,
                    std::vector<unsigned int> & global_hits_v,
                    unsigned int & global_hits_count) -> void
  {
    /* update swarm stats */
    auto const & seed_info = ampinfo_v[seed];

    ++current_swarm.size;
    current_swarm.maxgen = std::max(seed_info.generation, current_swarm.maxgen);
    const auto abundance = data.abundance(seed);
    current_swarm.abundance_sum += abundance;
    if (abundance == 1) {
      ++current_swarm.singletons;
    }
    current_swarm.sumlen += data.sequence_view(seed).length;

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

    for (auto offset = 0U; offset < link_count; ++offset)
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


  inline auto add_amp_to_swarm(unsigned int const amp,
                               std::vector<struct ampinfo_s> & ampinfo_v) -> void
  {
    /* add to swarm */
    ampinfo_v[current_swarm.tail].next = amp;
    current_swarm.tail = amp;
  }


  auto process_generation(unsigned int subseed,
                          Data const & data,
                          std::vector<struct ampinfo_s> & ampinfo_v,
                          std::vector<unsigned int> const & network_v,
                          std::vector<unsigned int> & global_hits_v) -> unsigned int
  {
    /* process all subseeds of this generation */
    auto global_hits_count = 0U;
    while (subseed != no_swarm)
      {
        process_seed(data, subseed, ampinfo_v, network_v, global_hits_v, global_hits_count);
        subseed = ampinfo_v[subseed].next;
      }

    /* sort all of this generation */
    std::sort(global_hits_v.begin(), global_hits_v.begin() + global_hits_count);

    /* add them to the swarm */
    for (auto i = 0U; i < global_hits_count; ++i) {
      add_amp_to_swarm(global_hits_v[i], ampinfo_v);
    }

    /* most abundant amplicon of next generation, or no_swarm if generation was empty */
    if (global_hits_count != 0U) {
      return global_hits_v[0];
    }
    return no_swarm;
  }


  auto ensure_swarm_capacity(unsigned int const swarmcount,
                             std::vector<struct swarminfo_s> & swarminfo_v) -> void
  {
    if (swarmcount >= swarminfo_v.size())
      {
        /* allocate memory for more swarms... */
        // note: capacity doubles, as usual
        // 1,024 times struct size (so at least 40,960 new bytes reserved)
        swarminfo_v.resize(swarminfo_v.size() + one_kilobyte);
      }
  }


  auto finalize_swarm_info(unsigned int const seed,
                           unsigned int const swarmcount,
                           std::vector<struct swarminfo_s> & swarminfo_v) -> void
  {
    auto & swarm_info = swarminfo_v[swarmcount];

    swarm_info.seed = seed;
    swarm_info.size = current_swarm.size;
    swarm_info.mass = current_swarm.abundance_sum;
    swarm_info.sumlen = current_swarm.sumlen;
    swarm_info.singletons = current_swarm.singletons;
    swarm_info.maxgen = current_swarm.maxgen;
    swarm_info.last = current_swarm.tail;
    swarm_info.attached = false;

    /* update overall stats */
    overall_stats.largest = std::max(current_swarm.size, overall_stats.largest);
    overall_stats.maxgen = std::max(current_swarm.maxgen, overall_stats.maxgen);
  }


  auto grow_swarm(unsigned int const seed,
                  unsigned int const swarmcount,
                  Data const & data,
                  std::vector<struct ampinfo_s> & ampinfo_v,
                  std::vector<struct swarminfo_s> & swarminfo_v,
                  std::vector<unsigned int> const & network_v,
                  std::vector<unsigned int> & global_hits_v) -> void
  {
    /* start a new swarm with a new initial seed */
    auto & seed_info = ampinfo_v[seed];
    seed_info.swarmid = swarmcount;
    seed_info.generation = 0;
    seed_info.parent = no_swarm;
    seed_info.next = no_swarm;

    /* initialize swarm stats and link up this initial seed in
       the list of swarms */
    current_swarm = Active_swarm_stats {};
    current_swarm.tail = seed;

    /* walk generations: each call processes the .next chain
       starting at subseed and returns the most abundant amplicon
       of the next generation (or no_swarm when exhausted). The
       first iteration starts from seed itself, whose .next is
       no_swarm, so it processes only the initial seed. */
    auto subseed = seed;
    while (subseed != no_swarm)
      {
        subseed = process_generation(subseed, data, ampinfo_v, network_v, global_hits_v);
      }

    ensure_swarm_capacity(swarmcount, swarminfo_v);
    finalize_swarm_info(seed, swarmcount, swarminfo_v);
  }


  auto write_network_file(const unsigned int number_of_networks,
                          struct Parameters const & parameters,
                          Data const & data,
                          std::vector<struct ampinfo_s> const & ampinfo_v,
                          std::vector<unsigned int> & network_v) -> void {
    // a network is a cluster with at least two sequences (no singletons)
    Progress progress("Dumping network:  ", number_of_networks, parameters);

    assert(ampinfo_v.size() == amplicons);
    auto counter = 0ULL;
    for (auto const& amplicon: ampinfo_v) {
      const auto link_start = amplicon.link_start;
      const auto link_count = amplicon.link_count;

      // amplicon indexes are already sorted by decreasing abundance
      // then by header in db.cc, so a natural ascending sort here
      // emits neighbours in that ranking order. Earlier dereplication
      // guarantees indexes are distinct.
      std::sort(network_v.begin() + link_start,
                network_v.begin() + link_start + link_count);

      // refactoring: std::vector<unsigned int> network_v(network_v.begin() + link_start, network_v.begin() + link_start + link_count);
      for (auto link = 0U; link < link_count; ++link)
        {
          const auto neighbour = network_v[link_start + link];
          data.fprint_id(parameters.network_file.get(), counter, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.network_file.get(), "\t");
          data.fprint_id(parameters.network_file.get(), neighbour, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.network_file.get(), "\n");
          progress.update();
        }
      ++counter;
    }
    progress.done();
  }


  auto write_swarms_default_format(struct Parameters const & parameters,
                                   Data const & data,
                                   std::vector<struct ampinfo_s> const & ampinfo_v,
                                   std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    static constexpr char sepchar {' '};
    Progress progress("Writing swarms:   ", swarminfo_v.size(), parameters);

    for (auto i = 0U; i < swarminfo_v.size(); ++i) {
      if (swarminfo_v[i].attached) {
        continue;
      }

      const auto seed = swarminfo_v[i].seed;
      for (auto amp_id = seed; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next) {
        if (amp_id != seed) {
          std::fputc(sepchar, parameters.outfile.get());
        }
        data.fprint_id(parameters.outfile.get(), amp_id,
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      }
      std::fputc('\n', parameters.outfile.get());
      progress.update(i + 1);
    }

    progress.done();
  }


  auto write_swarms_mothur_format(struct Parameters const & parameters,
                                  Data const & data,
                                  std::vector<struct ampinfo_s> const & ampinfo_v,
                                  std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    Progress progress("Writing swarms:   ", swarminfo_v.size(), parameters);

    std::fprintf(parameters.outfile.get(), "swarm_%" PRId64 "\t%" PRIu64,
                 parameters.opt_differences, overall_stats.swarmcount_adjusted);

    for (auto i = 0U; i < swarminfo_v.size(); ++i) {
      assert(not swarminfo_v[i].attached);
      if (swarminfo_v[i].attached) {
        continue;
      }

      const auto seed = swarminfo_v[i].seed;
      for (auto amp_id = seed; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next) {
        if (amp_id == seed) {
          std::fputc('\t', parameters.outfile.get());
        }
        else {
          std::fputc(',', parameters.outfile.get());
        }
        data.fprint_id(parameters.outfile.get(), amp_id,
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      }
      progress.update(i + 1);
    }

    std::fputc('\n', parameters.outfile.get());

    progress.done();
  }


  auto write_swarms_uclust_format(struct Parameters const & parameters,
                                  Data const & data,
                                  std::vector<struct ampinfo_s> const & ampinfo_v,
                                  std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    auto cluster_no = 0U;
    NwAligner aligner(data.longest_sequence(),
                      parameters.penalty_mismatch,
                      static_cast<unsigned long int>(parameters.penalty_gapopen),
                      static_cast<unsigned long int>(parameters.penalty_gapextend));

    Progress progress("Writing UCLUST:   ", swarminfo_v.size(), parameters);

    for (auto const & swarm_info : swarminfo_v) {
      if (swarm_info.attached) {
        continue;
      }

      const auto seed = swarm_info.seed;

      auto const & seed_info = ampinfo_v[seed];

      std::fprintf(parameters.uclustfile.get(), "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                   cluster_no,
                   swarm_info.size);
      data.fprint_id(parameters.uclustfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile.get(), "\t*\n");

      std::fprintf(parameters.uclustfile.get(), "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                   cluster_no,
                   data.sequence_view(seed).length);
      data.fprint_id(parameters.uclustfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile.get(), "\t*\n");

      for (auto amp_id = seed_info.next; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next)
        {
          auto const amp_seq = data.sequence_view(amp_id);
          auto const seed_seq = data.sequence_view(seed);  // refactoring: can be moved outside of this loop!

          auto const result = aligner.align(amp_seq, seed_seq);

          std::fprintf(parameters.uclustfile.get(),
                       "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                       cluster_no,
                       amp_seq.length,
                       result.percent_id,
                       result.differences > 0 ? result.cigar_string.data() : "=");

          data.fprint_id(parameters.uclustfile.get(), amp_id, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile.get(), "\t");
          data.fprint_id(parameters.uclustfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile.get(), "\n");
        }

      ++cluster_no;
      progress.update();
    }
    progress.done();
  }


  auto write_representative_sequences(struct Parameters const & parameters,
                                      Data const & data,
                                      std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    Progress progress("Writing seeds:    ", swarminfo_v.size(), parameters);

    std::vector<unsigned int> sorter(swarminfo_v.size());
    std::iota(sorter.begin(), sorter.end(), 0);

    auto compare_mass_and_headers = [&swarminfo_v, &data](unsigned int const lhs,
                                                          unsigned int const rhs) -> bool
    {
      const auto & swarm_x = swarminfo_v[lhs];
      const auto & swarm_y = swarminfo_v[rhs];

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
      // assert(data.header_view(swarm_x.seed) != data.header_view(swarm_y.seed)); // all headers are unique
      return data.header_view(swarm_x.seed) < data.header_view(swarm_y.seed);
    };

    std::sort(sorter.begin(), sorter.end(), compare_mass_and_headers);

    for (const auto index : sorter) {
      const auto & a_swarm = swarminfo_v[index];
      if (a_swarm.attached) {
        continue;
      }
      const auto seed = a_swarm.seed;
      const auto mass = a_swarm.mass;
      std::fprintf(parameters.seeds_file.get(), ">");
      data.fprint_id_with_new_abundance(parameters.seeds_file.get(), seed, mass,
                                   parameters.opt_usearch_abundance);
      std::fprintf(parameters.seeds_file.get(), "\n");
      data.fprintseq(parameters.seeds_file.get(), seed);
      progress.update();
    }

    progress.done();
  }


  auto write_structure_file(struct Parameters const & parameters,
                            Data const & data,
                            std::vector<struct ampinfo_s> const & ampinfo_v,
                            std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    auto cluster_no = 0U;

    Progress progress("Writing structure:", swarminfo_v.size(), parameters);

    for (auto swarmid = 0U; swarmid < swarminfo_v.size(); ++swarmid)
      {
        if (swarminfo_v[swarmid].attached) {
          continue;
        }
        const auto seed = swarminfo_v[swarmid].seed;

        auto const & seed_info = ampinfo_v[seed];

        for (auto amp_id = seed_info.next; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next)
          {
            const auto graft_parent = ampinfo_v[amp_id].graft_cand;
            if (graft_parent != no_swarm)
              {
                data.fprint_id_noabundance(parameters.internal_structure_file.get(),
                                      graft_parent, parameters.opt_usearch_abundance);
                std::fprintf(parameters.internal_structure_file.get(), "\t");
                data.fprint_id_noabundance(parameters.internal_structure_file.get(), amp_id, parameters.opt_usearch_abundance);
                std::fprintf(parameters.internal_structure_file.get(),
                             "\t%d\t%u\t%u\n",
                             2,
                             cluster_no + 1,
                             ampinfo_v[graft_parent].generation + 1);
              }

            const auto parent = ampinfo_v[amp_id].parent;
            if (parent != no_swarm)
              {
                data.fprint_id_noabundance(parameters.internal_structure_file.get(), parent, parameters.opt_usearch_abundance);
                std::fprintf(parameters.internal_structure_file.get(), "\t");
                data.fprint_id_noabundance(parameters.internal_structure_file.get(), amp_id, parameters.opt_usearch_abundance);
                std::fprintf(parameters.internal_structure_file.get(),
                             "\t%u\t%u\t%u\n",
                             1U,
                             cluster_no + 1,
                             ampinfo_v[amp_id].generation);
              }
          }

        ++cluster_no;
        progress.update(swarmid);
      }
    progress.done();
  }


  auto write_stats_file(struct Parameters const & parameters,
                        Data const & data,
                        std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    Progress progress("Writing stats:    ", swarminfo_v.size(), parameters);

    for (auto const & swarm_info : swarminfo_v) {
      assert(not swarm_info.attached);
      if (swarm_info.attached) {
        continue;
      }
      std::fprintf(parameters.statsfile.get(), "%u\t%" PRIu64 "\t", swarm_info.size, swarm_info.mass);
      data.fprint_id_noabundance(parameters.statsfile.get(), swarm_info.seed, parameters.opt_usearch_abundance);
      std::fprintf(parameters.statsfile.get(), "\t%" PRIu64 "\t%u\t%u\t%u\n",
                   data.abundance(swarm_info.seed),
                   swarm_info.singletons, swarm_info.maxgen, swarm_info.maxgen);
      progress.update();
    }
    progress.done();
  }


  auto output_results(struct Parameters const & parameters,
                      Data const & data,
                      std::vector<struct ampinfo_s> const & ampinfo_v,
                      std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    /* dump swarms */
    if (parameters.opt_mothur) {
      write_swarms_mothur_format(parameters, data, ampinfo_v, swarminfo_v);
    }
    else {
      write_swarms_default_format(parameters, data, ampinfo_v, swarminfo_v);
    }

    /* dump seeds in fasta format with sum of abundances */
    if (not parameters.opt_seeds.empty()) {
      write_representative_sequences(parameters, data, swarminfo_v);
    }

    /* output internal structure */
    if (not parameters.opt_internal_structure.empty()) {
      write_structure_file(parameters, data, ampinfo_v, swarminfo_v);
    }

    /* output swarms in uclust format */
    if (not parameters.opt_uclust_file.empty()) {
      write_swarms_uclust_format(parameters, data, ampinfo_v, swarminfo_v);
    }

    /* output statistics to file */
    if (not parameters.opt_statistics_file.empty()) {
      write_stats_file(parameters, data, swarminfo_v);
    }
  }


  auto count_cluster_stats(struct Parameters const & parameters,
                           std::vector<struct swarminfo_s> const & swarminfo_v) -> Cluster_stats
  {
    Cluster_stats stats;

    Progress progress_count("Counting amplicons in heavy and light swarms",
                            swarminfo_v.size(), parameters);

    for (auto const & swarm_info : swarminfo_v)
      {
        if (swarm_info.mass < static_cast<uint64_t>(parameters.opt_boundary))
          {
            stats.amplicons_in_small_clusters += swarm_info.size;
            stats.nucleotides_in_small_clusters += swarm_info.sumlen;
            ++stats.small_clusters;
          }
        progress_count.update();
      }
    progress_count.done();

    stats.amplicons_in_large_clusters = amplicons - stats.amplicons_in_small_clusters;
    stats.large_clusters = swarminfo_v.size() - stats.small_clusters;

    return stats;
  }


  auto compute_bloom_geometry(struct Parameters const & parameters,
                              uint64_t const nucleotides_in_small_clusters) -> Bloom_geometry
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

    auto bloom_length_in_bits = nucleotides_in_small_clusters * microvariants * bits;

    auto const memtotal = system_get_memtotal();
    auto const memused = system_get_memused();

    if (parameters.opt_ceiling != 0)
      {
        if (static_cast<uint64_t>(parameters.opt_ceiling) * one_megabyte < memused)
          {
            fatal("Memory ceiling for Bloom filter is too low.");
          }
        assert(memused < one_megabyte * static_cast<uint64_t>(parameters.opt_ceiling));
        const uint64_t memrest
          = (one_megabyte * static_cast<uint64_t>(parameters.opt_ceiling)) - memused;
        auto const new_bits = n_bits_in_a_byte * memrest / (microvariants * nucleotides_in_small_clusters);
        if (new_bits < bits)
          {
            if (new_bits < 2) {
              fatal("Insufficient memory remaining for Bloom filter.");
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

    if (memused + (bloom_length_in_bits / n_bits_in_a_byte) > memtotal)
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

    Bloom_geometry geom;
    geom.n_bytes = ((bloom_length_in_bits - 1) / n_bits_in_a_byte) + 1;
    geom.n_hash_functions = n_hash_functions;
    return geom;
  }


  auto run_light_pass(struct Parameters const & parameters,
                      Data const & data,
                      std::vector<struct ampinfo_s> const & ampinfo_v,
                      std::vector<struct swarminfo_s> const & swarminfo_v,
                      Hashtable & hash_table,
                      BloomFilter & bloom_a,
                      BloomFilter & bloom_f,
                      uint64_t const amplicons_in_small_clusters) -> void
  {
    Progress progress_light("Adding light swarm amplicons to Bloom filter",
                            amplicons_in_small_clusters, parameters);

    /* process amplicons in order from least to most abundant */
    /* but stop when all amplicons in small clusters are processed */

    struct Light_state light_state;
    light_state.amplicon_count = amplicons_in_small_clusters;
    light_state.amplicon = amplicons - 1;
    {
      auto const light_tr = utils::make_unique<ThreadRunner>(
          static_cast<std::size_t>(parameters.opt_threads),
          [&parameters, &data, &ampinfo_v, &swarminfo_v, &hash_table, &bloom_a, &bloom_f, &light_state, &progress_light](uint64_t nth_thread) -> void {
            mark_light_thread(parameters, data, ampinfo_v, swarminfo_v, hash_table, bloom_a, bloom_f, nth_thread, light_state, progress_light);
          });
      light_tr->run();
    }
    progress_light.done();

    std::fprintf(parameters.logfile,
                 "Generated %" PRIu64 " variants from light swarms\n",
                 light_state.variants);
  }


  auto run_heavy_pass(struct Parameters const & parameters,
                      Data const & data,
                      std::vector<struct ampinfo_s> & ampinfo_v,
                      std::vector<struct swarminfo_s> const & swarminfo_v,
                      Hashtable const & hash_table,
                      BloomFilter const & bloom_a,
                      BloomFilter const & bloom_f,
                      uint64_t const amplicons_in_large_clusters,
                      struct Graft_state & graft_state) -> void
  {
    Progress progress_heavy("Checking heavy swarm amplicons against Bloom filter",
                            amplicons_in_large_clusters, parameters);

    /* process amplicons in order from most to least abundant */
    /* but stop when all amplicons in large clusters are processed */

    struct Heavy_state heavy_state;
    heavy_state.amplicon_count = amplicons_in_large_clusters;
    {
      auto const heavy_tr = utils::make_unique<ThreadRunner>(
          static_cast<std::size_t>(parameters.opt_threads),
          [&parameters, &data, &ampinfo_v, &swarminfo_v, &hash_table, &bloom_a, &bloom_f, &heavy_state, &graft_state, &progress_heavy](uint64_t nth_thread) -> void {
            check_heavy_thread(parameters, data, ampinfo_v, swarminfo_v, hash_table, bloom_a, bloom_f, nth_thread, heavy_state, graft_state, progress_heavy);
          });
      heavy_tr->run();
    }
    progress_heavy.done();

    std::fprintf(parameters.logfile, "Heavy variants: %" PRIu64 "\n", heavy_state.variants);
    std::fprintf(parameters.logfile, "Got %" PRId64 " graft candidates\n", graft_state.candidates);
  }


  auto log_swarm_summary(struct Parameters const & parameters) -> void
  {
    std::fprintf(parameters.logfile, "\n");
    std::fprintf(parameters.logfile, "Number of swarms:  %" PRIu64 "\n", overall_stats.swarmcount_adjusted);
    std::fprintf(parameters.logfile, "Largest swarm:     %u\n", overall_stats.largest);
    std::fprintf(parameters.logfile, "Max generations:   %u\n", overall_stats.maxgen);
  }


  auto run_fastidious_pass(struct Parameters const & parameters,
                           Data const & data,
                           unsigned int const swarmcount,
                           std::vector<struct ampinfo_s> & ampinfo_v,
                           std::vector<struct swarminfo_s> & swarminfo_v) -> void
  {
    std::fprintf(parameters.logfile, "\n");
    std::fprintf(parameters.logfile, "Results before fastidious processing:\n");
    std::fprintf(parameters.logfile, "Number of swarms:  %u\n", swarmcount);
    std::fprintf(parameters.logfile, "Largest swarm:     %u\n", overall_stats.largest);
    std::fprintf(parameters.logfile, "\n");

    auto const stats = count_cluster_stats(parameters, swarminfo_v);
    auto const small_clusters = stats.small_clusters;
    auto const large_clusters = stats.large_clusters;
    auto const amplicons_in_small_clusters = stats.amplicons_in_small_clusters;
    auto const amplicons_in_large_clusters = stats.amplicons_in_large_clusters;
    auto const nucleotides_in_small_clusters = stats.nucleotides_in_small_clusters;

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
        auto const bloom_geom = compute_bloom_geometry(parameters, nucleotides_in_small_clusters);
        static constexpr unsigned int fastidious_pattern_shift {16};
        BloomFilter bloom_f(bloom_geom.n_bytes, fastidious_pattern_shift,
                            bloom_geom.n_hash_functions);


        /* Allocate a fresh per-amplicon hash table and Bloom filter
           for the fastidious phase; only light-cluster amplicons will
           be inserted. */
        Hashtable hash_table;
        const auto hashtablesize = hash_table.allocate(amplicons);
        BloomFilter bloom_a(hashtablesize, amplicon_pattern_shift,
                            amplicon_n_hash_functions);

        run_light_pass(parameters, data, ampinfo_v, swarminfo_v, hash_table, bloom_a, bloom_f, amplicons_in_small_clusters);

        struct Graft_state graft_state;
        run_heavy_pass(parameters, data, ampinfo_v, swarminfo_v, hash_table, bloom_a, bloom_f, amplicons_in_large_clusters, graft_state);

        auto const grafts = attach_candidates(parameters, amplicons, ampinfo_v, swarminfo_v);
        std::fprintf(parameters.logfile, "Made %u grafts\n", grafts);
        std::fprintf(parameters.logfile, "\n");
      }
  }
} // namespace


auto algo_d1_run(struct Parameters const & parameters,
                 Data const & data) -> void
{
  amplicons = data.sequence_count();

  std::vector<struct ampinfo_s> ampinfo_v(amplicons);

  std::vector<struct swarminfo_s> swarminfo_v(one_kilobyte);

  // max number of microvariants = 7 * len + 4
  static constexpr auto multiplier = 7U;
  static constexpr auto offset = 4U;
  const auto global_hits_alloc = (multiplier * data.longest_sequence()) + offset + 1;
  std::vector<unsigned int> global_hits_v(global_hits_alloc);
  global_hits_data = global_hits_v.data();


  /* for all amplicons, generate list of matching amplicons */
  struct Network_state network_state;
  network_state.network_v.resize(one_megabyte);

  /* d=1 hashtable and Bloom filter live in their own scope so their
     backing storage is released before run_fastidious_pass() allocates
     its own fresh pair. */
  {
    /* populate the d=1 hash table and Bloom filter with the amplicon
       hashes precomputed in db.cc */
    Hashtable hash_table;
    const auto hashtablesize = hash_table.allocate(amplicons);
    BloomFilter bloom_a(hashtablesize, amplicon_pattern_shift,
                        amplicon_n_hash_functions);

    Progress progress_hash("Building hashtable:", amplicons, parameters);

    for (auto k = 0U; k < amplicons; ++k)
      {
        hash_insert(data, hash_table, bloom_a, k);
        progress_hash.update(k);
      }

    progress_hash.done();


    Progress progress_network("Building network: ", amplicons, parameters);
    {
      auto const network_tr = utils::make_unique<ThreadRunner>(
          static_cast<std::size_t>(parameters.opt_threads),
          [&parameters, &data, &ampinfo_v, &hash_table, &bloom_a, &network_state, &progress_network](uint64_t nth_thread) -> void {
            network_thread(parameters, data, ampinfo_v, hash_table, bloom_a, nth_thread, network_state, progress_network);
          });
      network_tr->run();
    }

    progress_network.done();
  }


  /* dump network to file */
  if (not parameters.opt_network_file.empty()) {
    write_network_file(network_state.count, parameters, data, ampinfo_v, network_state.network_v);
  }


  /* for each non-swarmed amplicon look for subseeds ... */

  auto swarmcount = 0U;  // refactoring: find a way to know swarmcount in advance?
  Progress progress_cluster("Clustering:       ", amplicons, parameters);

  for (auto seed = 0U; seed < amplicons; ++seed)
    {
      if (ampinfo_v[seed].swarmid == no_swarm)
        {
          grow_swarm(seed, swarmcount, data, ampinfo_v, swarminfo_v,
                     network_state.network_v, global_hits_v);
          ++swarmcount;
        }
      progress_cluster.update(seed + 1);
    }
  progress_cluster.done();

  global_hits_data = nullptr;

  network_state.network_v.clear();
  network_state.network_v.shrink_to_fit();

  overall_stats.swarmcount_adjusted = swarmcount;

  /* fastidious */
  if (parameters.opt_fastidious) {
    run_fastidious_pass(parameters, data, swarmcount, ampinfo_v, swarminfo_v);
  }

  // refactoring: trim vectors (remove allocated unused elements)
  // could it be done before the fastidious phase?
  swarminfo_v.resize(swarmcount);  // swarminfo_v's capacity can be twice too much
  swarminfo_v.shrink_to_fit();

  output_results(parameters, data, ampinfo_v, swarminfo_v);

  log_swarm_summary(parameters);
}
