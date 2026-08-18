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

#include "algod1.hpp"
#include "swarm.hpp"
#include "db.hpp"
#include "utils/progress.hpp"
#include "utils/algod1_internal.hpp"
#include "utils/algod1_fastidious.hpp"
#include "utils/algod1_network.hpp"
#include "utils/algod1_output.hpp"
#include "utils/algod1_statistics.hpp"
#include "utils/span.hpp"  // Span<>, make_span()
#include <algorithm>  // std::sort(), std::max()
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <functional>  // std::reference_wrapper
#include <vector>


namespace {

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


  /* the hit buffer's initial size and its growth step, in entries: one
     unit for both, so that neither is derived from an unrelated quantity */
  constexpr std::size_t hit_chunk {4UL * one_kilobyte};


  /* One generation's hits, filled into storage that outlives them: the
     buffer is allocated once per run and reused by every generation of
     every cluster, so storage_.size() is the *allocated* length and never
     the used one. The used length lives here, next to the storage it
     applies to, instead of travelling beside it as a second parameter.

     Constructed once per generation, so starting empty is structural
     rather than something a caller has to remember to reset. */
  class Hit_list {
  public:
    explicit Hit_list(std::vector<unsigned int> & storage) noexcept
      : storage_(storage) {}

    /* room for 'additional' more hits, grown in whole chunks. used_ is a
       std::size_t, the same type as storage_.size(), so the sum below
       cannot wrap; as an unsigned int count it could, skipping the growth
       and leaving push() to write out of bounds. */
    auto reserve_for(std::size_t const additional) -> void {
      auto const required = used_ + additional;
      if (required <= storage().size()) { return; }
      auto enlarged = storage().size();
      while (required > enlarged) {
        enlarged += hit_chunk;
      }
      storage().resize(enlarged);
    }

    /* an indexed write, deliberately: push_back() here is the shape of the
       regression recorded in commit e517c04 (a candidate-list fill, 25-35 %
       at d=2) */
    auto push(unsigned int const amp) noexcept -> void {
      assert(used_ < storage().size());  // reserve_for() precedes the writes
      storage()[used_] = amp;
      ++used_;
    }

    /* the hits collected so far, mutable because the caller sorts them in
       place. Invalidated by any later push() that grows the storage. */
    auto filled() noexcept -> Span<unsigned int> {
      return make_span(storage()).first(used_);
    }

  private:
    /* a reference_wrapper rather than a plain 'std::vector<unsigned int> &':
       a reference member would delete the assignment operator, which is what
       cppcoreguidelines-avoid-const-or-ref-data-members reports. The borrow
       is unchanged -- the caller still owns the storage -- and storage()
       hands the vector back so the members above read as before. */
    std::reference_wrapper<std::vector<unsigned int>> storage_;
    std::size_t used_ {0};

    auto storage() noexcept -> std::vector<unsigned int> & { return storage_; }
  };


  auto process_seed(Data const & data,
                    unsigned int const seed,
                    std::vector<struct ampinfo_s> & ampinfo_v,
                    std::vector<unsigned int> const & network_v,
                    Hit_list & hits,
                    Active_swarm_stats & current_swarm) -> void
  {
    /* update swarm stats */
    auto const & seed_info = ampinfo_v[seed];

    ++current_swarm.size;
    current_swarm.maxgen = std::max(seed_info.generation, current_swarm.maxgen);
    auto const abundance = data.abundance(seed);
    current_swarm.abundance_sum += abundance;
    if (abundance == 1) {
      ++current_swarm.singletons;
    }
    current_swarm.sumlen += data.sequence_view(seed).length;

    auto const neighbours = neighbours_of(network_v, seed_info);
    hits.reserve_for(neighbours.size());

    for (auto const amp : neighbours)
      {
        if (ampinfo_v[amp].swarmid == no_swarm)
          {
            /* each amplicon is appended at most once in the whole run: the
               guard above is answered by the stamp below, and a stamped
               amplicon is never a hit again. The buffer therefore never
               needs more than one entry per amplicon. */
            hits.push(amp);

            /* update info */
            ampinfo_v[amp].swarmid = ampinfo_v[seed].swarmid;
            ampinfo_v[amp].generation = ampinfo_v[seed].generation + 1;
            ampinfo_v[amp].parent = seed;
          }
      }
  }


  inline auto add_amp_to_swarm(unsigned int const amp,
                               std::vector<struct ampinfo_s> & ampinfo_v,
                               Active_swarm_stats & current_swarm) -> void
  {
    /* add to swarm */
    ampinfo_v[current_swarm.tail].next = amp;
    current_swarm.tail = amp;
  }


  auto process_generation(unsigned int subseed,
                          Data const & data,
                          std::vector<struct ampinfo_s> & ampinfo_v,
                          std::vector<unsigned int> const & network_v,
                          std::vector<unsigned int> & global_hits_v,
                          Active_swarm_stats & current_swarm) -> unsigned int
  {
    /* process all subseeds of this generation */
    Hit_list hits {global_hits_v};
    while (subseed != no_swarm)
      {
        process_seed(data, subseed, ampinfo_v, network_v, hits, current_swarm);
        subseed = ampinfo_v[subseed].next;
      }

    /* sort all of this generation */
    auto const generation_hits = hits.filled();

    /* each amplicon is appended at most once in the whole run (see
       process_seed), and this cluster's own seed was stamped by grow_swarm
       without being appended, so one generation cannot hold as many hits as
       the database has amplicons. This is the bound that would let the
       buffer be pre-sized instead of grown; asserting it here costs nothing
       in release builds. */
    assert(generation_hits.size() < data.sequence_count());

    std::sort(generation_hits.begin(), generation_hits.end());

    /* add them to the swarm */
    for (auto const amp : generation_hits) {
      add_amp_to_swarm(amp, ampinfo_v, current_swarm);
    }

    /* most abundant amplicon of next generation, or no_swarm if generation was empty */
    if (not generation_hits.empty()) {
      return generation_hits.front();
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
                           std::vector<struct swarminfo_s> & swarminfo_v,
                           Active_swarm_stats const & current_swarm,
                           Overall_stats & overall_stats) -> void
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
                  std::vector<unsigned int> & global_hits_v,
                  Overall_stats & overall_stats) -> void
  {
    /* start a new swarm with a new initial seed */
    auto & seed_info = ampinfo_v[seed];
    seed_info.swarmid = swarmcount;
    seed_info.generation = 0;
    seed_info.parent = no_swarm;
    seed_info.next = no_swarm;

    /* initialize swarm stats and link up this initial seed in
       the list of swarms */
    Active_swarm_stats current_swarm {};
    current_swarm.tail = seed;

    /* walk generations: each call processes the .next chain
       starting at subseed and returns the most abundant amplicon
       of the next generation (or no_swarm when exhausted). The
       first iteration starts from seed itself, whose .next is
       no_swarm, so it processes only the initial seed. */
    auto subseed = seed;
    while (subseed != no_swarm)
      {
        subseed = process_generation(subseed, data, ampinfo_v, network_v, global_hits_v, current_swarm);
      }

    ensure_swarm_capacity(swarmcount, swarminfo_v);
    finalize_swarm_info(seed, swarmcount, swarminfo_v, current_swarm, overall_stats);
  }


  auto run_clustering(struct Parameters const & parameters,
                      Data const & data,
                      std::vector<struct ampinfo_s> & ampinfo_v,
                      std::vector<struct swarminfo_s> & swarminfo_v,
                      std::vector<unsigned int> const & network_v,
                      std::vector<unsigned int> & global_hits_v,
                      Overall_stats & overall_stats) -> unsigned int
  {
    /* for each non-swarmed amplicon look for subseeds ... */
    auto swarmcount = 0U;  // refactoring: find a way to know swarmcount in advance?
    auto const amplicons = data.sequence_count();
    Progress progress_cluster("Clustering:       ", amplicons, parameters);

    for (auto seed = 0U; seed < amplicons; ++seed)
      {
        if (ampinfo_v[seed].swarmid == no_swarm)
          {
            grow_swarm(seed, swarmcount, data, ampinfo_v, swarminfo_v,
                       network_v, global_hits_v, overall_stats);
            ++swarmcount;
          }
        progress_cluster.increment();
      }
    progress_cluster.done();
    return swarmcount;
  }

} // namespace


auto algo_d1_run(struct Parameters const & parameters,
                 Data const & data) -> void
{
  auto const amplicons = data.sequence_count();

  Overall_stats overall_stats {};

  std::vector<struct ampinfo_s> ampinfo_v(amplicons);

  std::vector<struct swarminfo_s> swarminfo_v(one_kilobyte);

  /* one growth chunk: the buffer is a scratch arena refilled once per
     generation and grown by Hit_list as needed, so its initial size is a
     starting point, not a bound. It used to be 7L+5, the microvariant upper
     bound for a single sequence, which sizes an unrelated quantity and fell
     below the growth step anyway. */
  std::vector<unsigned int> global_hits_v(hit_chunk);


  /* for all amplicons, generate list of matching amplicons */
  struct Network_state network_state;
  network_state.network_v.resize(one_megabyte);

  build_amplicon_network(parameters, data, ampinfo_v, network_state);


  /* dump network to file */
  if (not parameters.opt_network_file.empty()) {
    write_network_file(network_state.count, parameters, data, ampinfo_v, network_state.network_v);
  }


  auto const swarmcount = run_clustering(parameters, data, ampinfo_v, swarminfo_v,
                                         network_state.network_v, global_hits_v, overall_stats);

  network_state.network_v.clear();
  network_state.network_v.shrink_to_fit();

  overall_stats.swarmcount_adjusted = swarmcount;

  // swarminfo_v's size was grown in fixed-size chunks, so it can hold
  // default-constructed padding entries beyond swarmcount. Trim them
  // before the fastidious pass, which counts light swarms over the
  // whole vector: padding entries have mass 0 and would otherwise be
  // miscounted as light swarms (inflating the count and, with
  // --ceiling, dividing by zero in compute_bloom_geometry).
  swarminfo_v.resize(swarmcount);  // swarminfo_v's capacity can be twice too much
  swarminfo_v.shrink_to_fit();

  /* fastidious */
  if (parameters.opt_fastidious) {
    run_fastidious_pass(parameters, data, swarmcount, ampinfo_v, swarminfo_v, overall_stats);
  }

  output_results(parameters, data, ampinfo_v, swarminfo_v, overall_stats);

  log_swarm_summary(parameters, overall_stats);
}
