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
#include "db.h"
#include "utils/progress.h"
#include "utils/algod1_internal.h"
#include "utils/algod1_fastidious.h"
#include "utils/algod1_network.h"
#include "utils/algod1_output.h"
#include "utils/algod1_statistics.h"
#include <algorithm>  // std::sort(), std::max()
#include <cstdint>  // int64_t, uint64_t
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


  auto process_seed(Data const & data,
                    unsigned int const seed,
                    std::vector<struct ampinfo_s> & ampinfo_v,
                    std::vector<unsigned int> const & network_v,
                    std::vector<unsigned int> & global_hits_v,
                    unsigned int & global_hits_count,
                    Active_swarm_stats & current_swarm) -> void
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
    auto global_hits_count = 0U;
    while (subseed != no_swarm)
      {
        process_seed(data, subseed, ampinfo_v, network_v, global_hits_v, global_hits_count, current_swarm);
        subseed = ampinfo_v[subseed].next;
      }

    /* sort all of this generation */
    std::sort(global_hits_v.begin(), global_hits_v.begin() + global_hits_count);

    /* add them to the swarm */
    for (auto i = 0U; i < global_hits_count; ++i) {
      add_amp_to_swarm(global_hits_v[i], ampinfo_v, current_swarm);
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
    const auto amplicons = data.sequence_count();
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
  const auto amplicons = data.sequence_count();

  Overall_stats overall_stats {};

  std::vector<struct ampinfo_s> ampinfo_v(amplicons);

  std::vector<struct swarminfo_s> swarminfo_v(one_kilobyte);

  const auto global_hits_alloc = compute_microvariant_buffer_size(data.longest_sequence());
  std::vector<unsigned int> global_hits_v(global_hits_alloc);


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

  /* fastidious */
  if (parameters.opt_fastidious) {
    run_fastidious_pass(parameters, data, swarmcount, ampinfo_v, swarminfo_v, overall_stats);
  }

  // refactoring: trim vectors (remove allocated unused elements)
  // could it be done before the fastidious phase?
  swarminfo_v.resize(swarmcount);  // swarminfo_v's capacity can be twice too much
  swarminfo_v.shrink_to_fit();

  output_results(parameters, data, ampinfo_v, swarminfo_v, overall_stats);

  log_swarm_summary(parameters, overall_stats);
}
