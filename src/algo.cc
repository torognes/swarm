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

#include "swarm.h"
#include "db.h"
#include "utils/algo_internal.h"
#include "utils/qgram.h"
#include "utils/nw_aligner.h"
#include "scan.h"
#include "utils/make_unique.h"
#include "utils/progress.h"
#include "utils/search_data.h"
#include "utils/threads.h"  // ThreadRunner
#include <algorithm>  // std::min(), std::for_each
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc(), fflush
#include <cstdlib>  // qsort()
#include <iterator> // std::next
#include <limits>
#include <memory>  // unique pointer
#include <string>
#include <vector>


namespace {

  auto set_amplicon_ids(std::vector<struct ampliconinfo_s> & amplicons) -> void {
    // a simple id based on input order
    auto index = 0U;
    for (auto & amplicon : amplicons) {
      amplicon.ampliconid = index;
      ++index;
    }
  }


  auto set_bit_mode(struct Parameters const & parameters) -> int {
    static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
    static constexpr auto bit_mode_8 = 8;
    static constexpr auto bit_mode_16 = 16;

#ifdef __aarch64__
#if !defined(DEBUG) || !defined(COVERAGE)
    /* always use 16-bit version on aarch64 because it is faster */
    return bit_mode_16;
#endif
#endif

    // search in 16-bit mode when the number of allowed differences or
    // penalty values are high (8 bits are not enough to keep track of
    // the score)
    auto const diff_saturation
      = static_cast<uint64_t>(std::min(uint8_max / parameters.penalty_mismatch,
                                       uint8_max / (parameters.penalty_gapopen +
                                                    parameters.penalty_gapextend)));

    if (static_cast<uint64_t>(parameters.opt_differences) > diff_saturation) {
      return bit_mode_16;
    }
    return bit_mode_8;
  }


  auto collect_seeds(struct Parameters const & parameters,
                     Data const & data,
                     uint64_t const amplicons,
                     std::vector<struct ampliconinfo_s> const & amps_v) -> std::vector<struct swarminfo_t> {
    Progress progress("Collecting seeds:    ", amplicons, parameters);
    std::vector<struct swarminfo_t> seeds(amplicons);
    auto swarmcount = 0UL;
    uint64_t mass = 0;
    auto previous_id = amps_v[0].swarmid;
    auto seed = amps_v[0].ampliconid;
    mass += data.abundance(seed);
    for (auto i = 1ULL; i < amplicons; ++i) {
        auto const current_id = amps_v[i].swarmid;
        if (current_id != previous_id) {
            seeds[swarmcount].seed = seed;  // update previous
            seeds[swarmcount].mass = mass;
            ++swarmcount;
            mass = 0;
            seed = amps_v[i].ampliconid;
          }
        mass += data.abundance(amps_v[i].ampliconid);
        previous_id = current_id;
        progress.update(i);
      }
    seeds[swarmcount].seed = seed;
    seeds[swarmcount].mass = mass;
    ++swarmcount;

    // free some memory
    assert(swarmcount <= std::numeric_limits<long int>::max());
    seeds.erase(std::next(seeds.begin(), static_cast<long int>(swarmcount)), seeds.end());
    seeds.shrink_to_fit();

    return seeds;
  }


  auto sort_seeds(struct Parameters const & parameters,
                  Data const & data,
                  std::vector<struct swarminfo_t> & seeds) -> void {
    Progress const progress("Sorting seeds:    ", seeds.size(), parameters);

    auto compare_seeds = [&data](struct swarminfo_t const& lhs,
                                 struct swarminfo_t const& rhs) -> bool {
      // sort by decreasing mass...
      if (lhs.mass > rhs.mass) {
        return true;
      }
      if (lhs.mass < rhs.mass) {
        return false;
      }
      // ...then ties are sorted by label (alphabetical order)
      return data.header_view(lhs.seed) < data.header_view(rhs.seed);
    };

    std::sort(seeds.begin(), seeds.end(), compare_seeds);
    progress.done();
  }


  auto write_seeds(struct Parameters const & parameters,
                   Data const & data,
                   std::vector<struct swarminfo_t> const & seeds) -> void {
    Progress progress("Writing seeds:    ", seeds.size(), parameters);
    for (auto const& seed: seeds) {
      auto const swarm_mass = seed.mass;
      auto const swarm_seed = seed.seed;

      std::fprintf(parameters.seeds_file.get(), ">");
      data.fprint_id_with_new_abundance(parameters.seeds_file.get(), swarm_seed, swarm_mass, parameters.opt_usearch_abundance);
      std::fprintf(parameters.seeds_file.get(), "\n");
      data.fprintseq(parameters.seeds_file.get(), swarm_seed);
      progress.update();
    }
    progress.done();
  }


  auto find_correct_position_in_list(uint64_t const starting_position,
                                     uint64_t const target,
                                     uint64_t const seeded,
                                     struct ampliconinfo_s const & subseed,
                                     std::vector<struct ampliconinfo_s> const & amplicons) -> uint64_t {
    auto position = starting_position;  // == swarmed
    auto const targetampliconid = amplicons[target].ampliconid;

    while ((position > seeded) and
           (amplicons[position - 1].ampliconid > targetampliconid) and
           (amplicons[position - 1].generation > subseed.generation)) {
      --position;
    }
    return position;
  }


  auto move_target_to_first_unswarmed_position(uint64_t const position,
                                               uint64_t const target,
                                               std::vector<struct ampliconinfo_s> & amplicons) -> void {
    /* move the 'target' to the position ('swarmed')
       of the first unswarmed amplicon in the pool,
       then move the target further into the swarmed
       but unseeded part of the list, so that the
       swarmed amplicons are ordered by id */

    // position       target
    // |              |
    // A    B    C    D  (initial situation)
    //                   (rotate to the right)
    // D    A    B    C  (final situation)

    if (target <= position) { return; }
    // target > position
    // assert(target < amplicons.size());
    auto const temp = amplicons[target];  // refactoring: static?
    for (auto i = target; i > position; --i) {
      amplicons[i] = amplicons[i - 1];
    }
    amplicons[position] = temp;

    // Do not refactor with std::rotate: both reverse-iterator (2x slower,
    // reverted in commit 8662a4d) and forward-iterator (1.2x slower on
    // d=2 18SV9) variants regressed performance versus the raw loop.
  }


  auto write_representative_sequences(uint64_t const amplicons,
                                      struct Parameters const & parameters,
                                      Data const & data,
                                      std::vector<struct ampliconinfo_s> const & amps_v) -> void {
    auto seeds = collect_seeds(parameters, data, amplicons, amps_v);
    sort_seeds(parameters, data, seeds);
    write_seeds(parameters, data, seeds);
  }


  auto write_swarms_default_format(uint64_t const amplicons,
                                   struct Parameters const & parameters,
                                   Data const & data,
                                   std::vector<struct ampliconinfo_s> const & amps_v) -> void {
    /* native swarm output */
    static constexpr char sepchar {' '};  /* usually a space */
    static constexpr char sep_swarms {'\n'};

    data.fprint_id(parameters.outfile.get(), amps_v[0].ampliconid,
              parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    int64_t previous_id = amps_v[0].swarmid;

    for (auto i = 1ULL; i < amplicons; ++i) {
        int64_t const current_id = amps_v[i].swarmid;
        if (current_id == previous_id) {
          std::fputc(sepchar, parameters.outfile.get());
        }
        else {
          std::fputc(sep_swarms, parameters.outfile.get());
        }
        data.fprint_id(parameters.outfile.get(), amps_v[i].ampliconid,
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        previous_id = current_id;
      }
    std::fputc('\n', parameters.outfile.get());
  }


  auto write_swarms_mothur_format(uint64_t const amplicons,
                                  unsigned int const swarmid,
                                  struct Parameters const & parameters,
                                  Data const & data,
                                  std::vector<struct ampliconinfo_s> const & amps_v) -> void {
    /* mothur list file output */
    static constexpr char sep_amplicons {','};
    static constexpr char sep_swarms {'\t'};

    std::fprintf(parameters.outfile.get(), "swarm_%" PRId64 "\t%u\t", parameters.opt_differences, swarmid);

    data.fprint_id(parameters.outfile.get(), amps_v[0].ampliconid,
              parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    int64_t previous_id = amps_v[0].swarmid;

    for (auto i = 1ULL; i < amplicons; ++i) {
        int64_t const current_id = amps_v[i].swarmid;
        if (current_id == previous_id) {
          std::fputc(sep_amplicons, parameters.outfile.get());
        }
        else {
          std::fputc(sep_swarms, parameters.outfile.get());
        }
        data.fprint_id(parameters.outfile.get(), amps_v[i].ampliconid,
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        previous_id = current_id;
      }

    std::fputc('\n', parameters.outfile.get());
  }


  auto write_internal_structure_line(uint64_t const parent_id,
                                     uint64_t const child_id,
                                     uint64_t const diff,
                                     unsigned int const swarmid,
                                     unsigned int const generation,
                                     struct Parameters const & parameters,
                                     Data const & data) -> void {
    data.fprint_id_noabundance(parameters.internal_structure_file.get(),
                               parent_id, parameters.opt_usearch_abundance);
    std::fprintf(parameters.internal_structure_file.get(), "\t");
    data.fprint_id_noabundance(parameters.internal_structure_file.get(),
                               child_id, parameters.opt_usearch_abundance);
    std::fprintf(parameters.internal_structure_file.get(), "\t%" PRIu64, diff);
    std::fprintf(parameters.internal_structure_file.get(),
                 "\t%u\t%u\n",
                 swarmid, generation);
  }


  auto build_remaining_amplicons_list(uint64_t const swarmed,
                                      uint64_t const seed_abundance,
                                      struct Parameters const & parameters,
                                      Data const & data,
                                      std::vector<struct ampliconinfo_s> const & amps_v,
                                      Cluster_workspace & workspace) -> uint64_t {
    workspace.qgramamps_v.clear();
    std::for_each(std::next(amps_v.cbegin(), static_cast<long int>(swarmed)), amps_v.cend(),
        [&parameters, &data, seed_abundance, &workspace](
            struct ampliconinfo_s const & amplicon) -> void {
          auto const ampliconid = amplicon.ampliconid;
          if ((parameters.opt_no_cluster_breaking) or
              (data.abundance(ampliconid) <= seed_abundance)) {
            workspace.qgramamps_v.push_back(ampliconid);
          }
        });
    return workspace.qgramamps_v.size();
  }


  auto build_subseed_candidate_list(uint64_t const swarmed,
                                    uint64_t const amplicons,
                                    struct ampliconinfo_s const & subseed,
                                    struct Parameters const & parameters,
                                    Data const & data,
                                    std::vector<struct ampliconinfo_s> const & amps_v,
                                    Cluster_workspace & workspace) -> uint64_t {
    auto const subseed_abundance = data.abundance(subseed.ampliconid);
    workspace.qgramamps_v.clear();
    workspace.qgramindices_v.clear();
    for (auto i = swarmed; i < amplicons; ++i) {
      uint64_t const targetampliconid = amps_v[i].ampliconid;
      if ((amps_v[i].diffestimate <=
           subseed.radius + parameters.opt_differences) and
          ((parameters.opt_no_cluster_breaking) or
           (data.abundance(targetampliconid)
            <= subseed_abundance))) {
        workspace.qgramamps_v.push_back(targetampliconid);
        workspace.qgramindices_v.push_back(i);
      }
    }
    return workspace.qgramamps_v.size();
  }


  auto include_amplicon_in_cluster(uint64_t const position,
                                   uint64_t const diff,
                                   unsigned int const swarmid,
                                   struct ampliconinfo_s const & parent,
                                   std::vector<struct ampliconinfo_s> & amps_v,
                                   std::vector<uint64_t> & hits,
                                   Cluster_state & state,
                                   struct Parameters const & parameters,
                                   Data const & data) -> void {
    amps_v[position].swarmid = swarmid;
    assert(parent.generation <= std::numeric_limits<unsigned int>::max() - 1);
    amps_v[position].generation = parent.generation + 1;
    state.maxgen = std::max<uint64_t>(state.maxgen, amps_v[position].generation);
    assert(parent.radius <= std::numeric_limits<unsigned int>::max() - diff);
    amps_v[position].radius = static_cast<unsigned int>(parent.radius + diff);
    state.maxradius = std::max<uint64_t>(amps_v[position].radius, state.maxradius);

    auto const poolampliconid = amps_v[position].ampliconid;
    hits[state.hitcount] = poolampliconid;
    ++state.hitcount;

    if (not parameters.opt_internal_structure.empty()) {
      write_internal_structure_line(parent.ampliconid, poolampliconid, diff,
                                    swarmid, amps_v[position].generation,
                                    parameters, data);
    }

    auto const abundance = data.abundance(poolampliconid);
    state.amplicons_copies += abundance;
    if (abundance == 1) {
      ++state.singletons;
    }
    ++state.swarmsize;
  }


  auto write_stats_line(uint64_t const swarmsize,
                        uint64_t const amplicons_copies,
                        uint64_t const singletons,
                        uint64_t const maxgen,
                        uint64_t const maxradius,
                        uint64_t const seedampliconid,
                        struct Parameters const & parameters,
                        Data const & data) -> void {
    auto const abundance = data.abundance(seedampliconid);

    std::fprintf(parameters.statsfile.get(), "%" PRIu64 "\t%" PRIu64 "\t",
                 swarmsize, amplicons_copies);
    data.fprint_id_noabundance(parameters.statsfile.get(), seedampliconid, parameters.opt_usearch_abundance);
    std::fprintf(parameters.statsfile.get(),
                 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                 abundance, singletons, maxgen, maxradius);
  }


  auto write_uclust_cluster(unsigned int const swarmid,
                            uint64_t const swarmsize,
                            uint64_t const seedampliconid,
                            uint64_t const hitcount,
                            std::vector<uint64_t> const & hits,
                            NwAligner & aligner,
                            struct Parameters const & parameters,
                            Data const & data) -> void {
    std::fprintf(parameters.uclustfile.get(), "C\t%u\t%" PRIu64 "\t*\t*\t*\t*\t*\t",
            swarmid - 1, swarmsize);
    data.fprint_id(parameters.uclustfile.get(), seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    std::fprintf(parameters.uclustfile.get(), "\t*\n");

    std::fprintf(parameters.uclustfile.get(), "S\t%u\t%u\t*\t*\t*\t*\t*\t",
            swarmid - 1, data.sequence_view(seedampliconid).length);
    data.fprint_id(parameters.uclustfile.get(), seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    std::fprintf(parameters.uclustfile.get(), "\t*\n");
    std::fflush(parameters.uclustfile.get());

    for (auto i = 1ULL; i < hitcount; ++i) {
      auto const hit = hits[i];
      auto const hit_seq = data.sequence_view(hit);
      auto const seed_seq = data.sequence_view(seedampliconid);

      auto const result = aligner.align(hit_seq, seed_seq);

      std::fprintf(parameters.uclustfile.get(), "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                   swarmid - 1, hit_seq.length, result.percent_id,
                   result.differences > 0 ? result.cigar_string : "=");

      data.fprint_id(parameters.uclustfile.get(), hit, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile.get(), "\t");
      data.fprint_id(parameters.uclustfile.get(), seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile.get(), "\n");
      std::fflush(parameters.uclustfile.get());
    }
  }


  auto finalize_algo_run(uint64_t const amplicons,
                         unsigned int const swarmid,
                         uint64_t const largestswarm,
                         uint64_t const maxgenerations,
                         struct Parameters const & parameters,
                         Data const & data,
                         std::vector<struct ampliconinfo_s> const & amps_v,
                         struct Search_state & search_state) -> void {
    /* output swarms */
    if (amplicons != 0) {
      if (parameters.opt_mothur) {
        write_swarms_mothur_format(amplicons, swarmid, parameters, data, amps_v);
      }
      else {
        write_swarms_default_format(amplicons, parameters, data, amps_v);
      }
    }


    /* dump seeds in fasta format with sum of abundances */
    if ((not parameters.opt_seeds.empty()) and (amplicons != 0)) {
      write_representative_sequences(amplicons, parameters, data, amps_v);
    }

    std::fprintf(parameters.logfile, "\n");

    std::fprintf(parameters.logfile, "Number of swarms:  %u\n", swarmid);

    std::fprintf(parameters.logfile, "Largest swarm:     %" PRIu64 "\n", largestswarm);

    std::fprintf(parameters.logfile, "Max generations:   %" PRIu64 "\n", maxgenerations);

    search_end(search_state);
  }


  auto start_new_cluster(Pool_cursor & cursor,
                         unsigned int const swarmid,
                         std::vector<struct ampliconinfo_s> & amps_v,
                         Cluster_state & state,
                         Cluster_workspace & workspace,
                         Data const & data) -> uint64_t {
    /* process each initial seed */
    workspace.qgramamps_v.clear();

    auto const seedindex = cursor.seeded;
    ++cursor.seeded;

    amps_v[seedindex].swarmid = swarmid;

    uint64_t const seedampliconid = amps_v[seedindex].ampliconid;
    workspace.hits[state.hitcount] = seedampliconid;
    ++state.hitcount;

    auto const abundance = data.abundance(seedampliconid);
    state.amplicons_copies += abundance;
    if (abundance == 1) {
      ++state.singletons;
    }

    ++cursor.swarmed;

    return seedampliconid;
  }


  auto seed_first_generation(Search_context const & ctx,
                             Pool_cursor & cursor,
                             unsigned int const swarmid,
                             uint64_t const seedindex,
                             std::vector<struct ampliconinfo_s> & amps_v,
                             Cluster_workspace & workspace,
                             Cluster_state & state) -> void {
    /* find diff estimates between seed and each amplicon in pool */
    uint64_t const seedampliconid = amps_v[seedindex].ampliconid;
    auto const seed_abundance = ctx.data.abundance(seedampliconid);

    uint64_t const listlen = build_remaining_amplicons_list(cursor.swarmed, seed_abundance,
                                                            ctx.parameters, ctx.data,
                                                            amps_v, workspace);

    ctx.qgram_differ.fast(seedampliconid, workspace.qgramamps_v, workspace.qgramdiffs_v);

    uint64_t targetcount = 0;
    for (auto i = 0ULL; i < listlen; ++i) {
      auto const poolampliconid = workspace.qgramamps_v[i];
      auto const diff = workspace.qgramdiffs_v[i];
      assert(diff <= std::numeric_limits<unsigned int>::max());
      amps_v[cursor.swarmed + i].diffestimate = static_cast<unsigned int>(diff);
      if (diff <= static_cast<uint64_t>(ctx.parameters.opt_differences)) {
        workspace.targetindices[targetcount] = cursor.swarmed + i;
        workspace.targetampliconids[targetcount] = poolampliconid;
        ++targetcount;
      }
    }

    if (targetcount == 0) { return; }

    search_do(ctx.parameters, ctx.data, ctx.search_state, seedampliconid, targetcount, workspace.targetampliconids.data(),
              workspace.scores_v.data(), workspace.diffs_v.data(), workspace.alignlengths.data(), ctx.bits, ctx.search_threads);

    for (auto target_id = 0ULL; target_id < targetcount; ++target_id) {
      auto const diff = workspace.diffs_v[target_id];

      if (diff > static_cast<uint64_t>(ctx.parameters.opt_differences)) { continue; }
      auto const target = workspace.targetindices[target_id];

      /* move the 'target' to the position ('swarmed')
         of the first unswarmed amplicon in the pool */
      move_target_to_first_unswarmed_position(cursor.swarmed, target, amps_v);

      include_amplicon_in_cluster(cursor.swarmed, diff, swarmid,
                                  amps_v[seedindex], amps_v, workspace.hits,
                                  state, ctx.parameters, ctx.data);
      ++cursor.swarmed;
    }
  }


  auto grow_cluster_from_subseeds(Search_context const & ctx,
                                  Pool_cursor & cursor,
                                  unsigned int const swarmid,
                                  uint64_t const amplicons,
                                  std::vector<struct ampliconinfo_s> & amps_v,
                                  Cluster_workspace & workspace,
                                  Cluster_state & state) -> void {
    while (cursor.seeded < cursor.swarmed) {

      /* process each subseed */

      auto const & subseed = amps_v[cursor.seeded];

      ++cursor.seeded;

      uint64_t targetcount = 0;

      auto const subseedlistlen = build_subseed_candidate_list(cursor.swarmed, amplicons,
                                                               subseed, ctx.parameters,
                                                               ctx.data, amps_v, workspace);

      ctx.qgram_differ.fast(subseed.ampliconid,
                            workspace.qgramamps_v, workspace.qgramdiffs_v);

      for (auto i = 0ULL; i < subseedlistlen; ++i) {
        if (workspace.qgramdiffs_v[i] <= static_cast<uint64_t>(ctx.parameters.opt_differences)) {
          workspace.targetindices[targetcount] = workspace.qgramindices_v[i];
          workspace.targetampliconids[targetcount] = workspace.qgramamps_v[i];
          ++targetcount;
        }
      }

      if (targetcount == 0) { continue; }

      search_do(ctx.parameters, ctx.data, ctx.search_state, subseed.ampliconid, targetcount, workspace.targetampliconids.data(),
                workspace.scores_v.data(), workspace.diffs_v.data(), workspace.alignlengths.data(), ctx.bits, ctx.search_threads);

      for (auto target_id = 0ULL; target_id < targetcount; ++target_id) {
        auto const diff = workspace.diffs_v[target_id];

        if (diff > static_cast<uint64_t>(ctx.parameters.opt_differences)) { continue; }
        auto const target = workspace.targetindices[target_id];

        /* find correct position in list */

        auto const pos = find_correct_position_in_list(cursor.swarmed, target, cursor.seeded,
                                                       subseed, amps_v);

        move_target_to_first_unswarmed_position(pos, target, amps_v);

        include_amplicon_in_cluster(pos, diff, swarmid, subseed,
                                    amps_v, workspace.hits, state,
                                    ctx.parameters, ctx.data);
        ++cursor.swarmed;
      }
    }
  }


  auto write_cluster_outputs(unsigned int const swarmid,
                             uint64_t const seedampliconid,
                             Cluster_state const & state,
                             std::vector<uint64_t> const & hits,
                             NwAligner * const aligner,
                             struct Parameters const & parameters,
                             Data const & data) -> void {
    if (parameters.uclustfile.get() != nullptr) {
      write_uclust_cluster(swarmid, state.swarmsize, seedampliconid, state.hitcount, hits,
                           *aligner, parameters, data);
    }

    if (parameters.statsfile.get() != nullptr) {
      write_stats_line(state.swarmsize, state.amplicons_copies, state.singletons,
                       state.maxgen, state.maxradius, seedampliconid, parameters, data);
    }
  }
} // namespace


auto algo_run(struct Parameters const & parameters,
              Data const & data) -> void {
  std::vector<struct Search_data> search_data_v(static_cast<uint64_t>(parameters.opt_threads));
  struct Search_state search_state;
  search_begin(parameters, data, search_state, search_data_v);
  /* start threads */
  auto const search_threads = utils::make_unique<ThreadRunner>(
      static_cast<std::size_t>(parameters.opt_threads),
      [&parameters, &data, &search_state](uint64_t thread_id) -> void {
        search_worker_core(parameters, data, thread_id, search_state);
      });

  uint64_t largestswarm {0};
  uint64_t maxgenerations {0};

  auto const amplicons = data.sequence_count();
  uint64_t const longestamplicon = data.longest_sequence();

  // RAII: builds the qgram store, ThreadRunner is destroyed (workers
  // joined) at end of algo_run scope.
  QgramDiffer qgram_differ(parameters, data);

  std::vector<struct ampliconinfo_s> amps_v(amplicons);
  Cluster_workspace workspace(amplicons);

  // NwAligner is only needed when UCLUST output is requested; its
  // scratch buffers grow with longestamplicon^2, so allocate lazily.
  // The score matrix is built inside the constructor as well.
  std::unique_ptr<NwAligner> aligner;
  if (parameters.uclustfile.get() != nullptr) {
    aligner = utils::make_unique<NwAligner>(
        longestamplicon,
        parameters.penalty_mismatch,
        static_cast<unsigned long int>(parameters.penalty_gapopen),
        static_cast<unsigned long int>(parameters.penalty_gapextend));
  }

  set_amplicon_ids(amps_v);
  auto const bits = set_bit_mode(parameters);

  Search_context const ctx {parameters, data, qgram_differ,
                            search_state, search_threads.get(), bits};

  Pool_cursor cursor;

  auto swarmid = 0U;

  Progress progress("Clustering:       ", amplicons, parameters);
  while (cursor.seeded < amplicons) {

      ++swarmid;
      Cluster_state state;
      uint64_t const seedindex = cursor.seeded;

      uint64_t const seedampliconid = start_new_cluster(cursor, swarmid, amps_v,
                                                        state, workspace, data);

      seed_first_generation(ctx, cursor, swarmid, seedindex,
                            amps_v, workspace, state);

      grow_cluster_from_subseeds(ctx, cursor, swarmid, amplicons,
                                 amps_v, workspace, state);

      largestswarm = std::max(state.swarmsize, largestswarm);
      maxgenerations = std::max(state.maxgen, maxgenerations);

      write_cluster_outputs(swarmid, seedampliconid, state, workspace.hits,
                            aligner.get(), parameters, data);

      progress.update(cursor.seeded);
  }
  progress.done();

  finalize_algo_run(amplicons, swarmid, largestswarm, maxgenerations,
                    parameters, data, amps_v, search_state);
}

