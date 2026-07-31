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

#include "algod1_fastidious.hpp"
#include "../swarm.hpp"
#include "../db.hpp"
#include "variants.hpp"
#include "algod1_internal.hpp"
#include "algod1_statistics.hpp"
#include "bloom.hpp"
#include "hashtable.hpp"
#include "make_unique.hpp"
#include "nt_codec.hpp"
#include "progress.hpp"
#include "threads.hpp"
#include "view.hpp"
#include <algorithm>  // std::sort(), std::max(), std::count_if()
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fprintf(), fputc(), fputs()
#include <mutex>  // std::lock_guard, std::unique_lock
#include <vector>


namespace {

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

  /******************** FASTIDIOUS START ********************/


  auto attach(unsigned int seed, unsigned int amp,
              std::vector<struct ampinfo_s> & ampinfo_v,
              std::vector<struct swarminfo_s> & swarminfo_v,
              Overall_stats & overall_stats) -> void
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


  auto count_pairs(std::vector<struct ampinfo_s> const & ampinfo_v) -> unsigned int {
    return static_cast<unsigned int>(
      std::count_if(ampinfo_v.cbegin(), ampinfo_v.cend(),
                    [](struct ampinfo_s const & info) -> bool {
                      return info.graft_cand != no_swarm;
                    }));
  }


  auto attach_candidates(struct Parameters const & parameters,
                         unsigned int amplicon_count,
                         std::vector<struct ampinfo_s> & ampinfo_v,
                         std::vector<struct swarminfo_s> & swarminfo_v,
                         Overall_stats & overall_stats) -> unsigned int
  {
    auto const pair_count = count_pairs(ampinfo_v);

    Progress progress("Grafting light swarms on heavy swarms", pair_count, parameters);

    /* allocate memory */
    std::vector<struct graft_cand> graft_array(pair_count);

    /* fill in */
    assert(ampinfo_v.size() == amplicon_count);
    auto ticker = 0U;
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
          attach(parent, child, ampinfo_v, swarminfo_v, overall_stats);
          ++grafts;
        }
      progress.increment();
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
    for (auto const & var : make_view(variant_list).first(variant_count)) {
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
            auto const var_seq = Sequence{make_view(varseq).first(nt_bytelength(varlen)), varlen};
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
                          struct Heavy_state & heavy_state,
                          struct Graft_state & graft_state,
                          Progress & progress) -> void
  {
    static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
    static constexpr auto offset = 4U;
    static constexpr auto nt_per_uint64 = 32U;  // 32 nucleotides can fit in a uint64

    std::vector<struct var_s> variant_list((multiplier * data.longest_sequence()) + offset);
    std::vector<struct var_s> variant_list2((multiplier * (data.longest_sequence() + 1)) + offset);

    const std::size_t size =
      sizeof(uint64_t) * ((data.longest_sequence() + 2 + nt_per_uint64 - 1) / nt_per_uint64);
    std::vector<char> buffer1(size);
    const auto amplicons = data.sequence_count();
    std::unique_lock<std::mutex> lock(heavy_state.mutex);
    while ((heavy_state.amplicon < amplicons) and
           (heavy_state.progress < heavy_state.amplicon_count))
      {
        auto const heavy_amplicon_id = heavy_state.amplicon;
        ++heavy_state.amplicon;
        auto const & target_amplicon = ampinfo_v[heavy_amplicon_id];
        auto const & target_swarm = swarminfo_v[target_amplicon.swarmid];
        if (target_swarm.mass >= parameters.opt_boundary)
          {
            ++heavy_state.progress;
            progress.update(heavy_state.progress);
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
                         struct Light_state & state,
                         Progress & progress) -> void
  {
    static constexpr auto multiplier = 7U;  // max number of microvariants = 7 * len + 4
    static constexpr auto offset = 4U;

    std::vector<struct var_s> variant_list((multiplier * data.longest_sequence()) + offset);

    std::unique_lock<std::mutex> lock(state.mutex);
    while (state.progress < state.amplicon_count)
      {
        const auto light_amplicon_id = state.amplicon;
        // Invariant: amplicon_count equals the number of light-swarm
        // amplicons in [0, sequence_count), so the loop stops before this
        // unsigned cursor underflows. Assert it so a future change to the
        // counting in count_cluster_stats() cannot silently become an
        // out-of-bounds read here.
        assert(light_amplicon_id < ampinfo_v.size());
        --state.amplicon;
        auto const & target_amplicon = ampinfo_v[light_amplicon_id];
        auto const & target_swarm = swarminfo_v[target_amplicon.swarmid];
        if (target_swarm.mass < parameters.opt_boundary)
          {
            ++state.progress;
            progress.update(state.progress);
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
    light_state.amplicon = data.sequence_count() - 1;
    {
      auto const light_tr = utils::make_unique<ThreadRunner>(
          parameters.opt_threads,
          [&parameters, &data, &ampinfo_v, &swarminfo_v, &hash_table, &bloom_a, &bloom_f, &light_state, &progress_light](uint64_t /*nth_thread*/) -> void {
            mark_light_thread(parameters, data, ampinfo_v, swarminfo_v, hash_table, bloom_a, bloom_f, light_state, progress_light);
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
          parameters.opt_threads,
          [&parameters, &data, &ampinfo_v, &swarminfo_v, &hash_table, &bloom_a, &bloom_f, &heavy_state, &graft_state, &progress_heavy](uint64_t /*nth_thread*/) -> void {
            check_heavy_thread(parameters, data, ampinfo_v, swarminfo_v, hash_table, bloom_a, bloom_f, heavy_state, graft_state, progress_heavy);
          });
      heavy_tr->run();
    }
    progress_heavy.done();

    std::fprintf(parameters.logfile, "Heavy variants: %" PRIu64 "\n", heavy_state.variants);
    std::fprintf(parameters.logfile, "Got %" PRId64 " graft candidates\n", graft_state.candidates);
  }

} // namespace


auto run_fastidious_pass(struct Parameters const & parameters,
                         Data const & data,
                         unsigned int const swarmcount,
                         std::vector<struct ampinfo_s> & ampinfo_v,
                         std::vector<struct swarminfo_s> & swarminfo_v,
                         Overall_stats & overall_stats) -> void
{
  const auto amplicons = data.sequence_count();

  static_cast<void>(std::fputc('\n', parameters.logfile));
  static_cast<void>(std::fputs("Results before fastidious processing:\n", parameters.logfile));
  std::fprintf(parameters.logfile, "Number of swarms:  %u\n", swarmcount);
  std::fprintf(parameters.logfile, "Largest swarm:     %u\n", overall_stats.largest);
  static_cast<void>(std::fputc('\n', parameters.logfile));

  auto const stats = count_cluster_stats(parameters, amplicons, swarminfo_v);
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

      auto const grafts = attach_candidates(parameters, amplicons, ampinfo_v, swarminfo_v, overall_stats);
      std::fprintf(parameters.logfile, "Made %u grafts\n", grafts);
      static_cast<void>(std::fputc('\n', parameters.logfile));
    }
}
