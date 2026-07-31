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

#include "algo.hpp"
#include "swarm.hpp"
#include "db.hpp"
#include "utils/algo_internal.hpp"
#include "utils/algo_output.hpp"
#include "utils/qgram.hpp"
#include "utils/nw_aligner.hpp"
#include "utils/scanner.hpp"
#include "utils/memory_budget.hpp"
#include "utils/make_unique.hpp"
#include "utils/progress.hpp"
#include <algorithm>  // std::min(), std::for_each
#include <cassert>
#include <cinttypes>  // macro PRIu64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fprintf(), fputc(), fflush
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


  auto set_bit_mode(struct Parameters const & parameters) -> Bit_mode {
    static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();

#ifdef __aarch64__
#if !defined(DEBUG) && !defined(COVERAGE)
    /* always use 16-bit version on aarch64 because it is faster */
    return Bit_mode::bits_16;
#endif
#endif

    // search in 16-bit mode when the number of allowed differences or
    // penalty values are high (8 bits are not enough to keep track of
    // the score)
    auto const diff_saturation
      = static_cast<uint64_t>(std::min(uint8_max / parameters.penalty_mismatch,
                                       uint8_max / (parameters.penalty_gapopen +
                                                    parameters.penalty_gapextend)));

    if (parameters.opt_differences > diff_saturation) {
      return Bit_mode::bits_16;
    }
    return Bit_mode::bits_8;
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
    auto const temp = amplicons[target];
    for (auto i = target; i > position; --i) {
      amplicons[i] = amplicons[i - 1];
    }
    amplicons[position] = temp;

    // Do not refactor with std::rotate: both reverse-iterator (2x slower,
    // reverted in commit 8662a4d) and forward-iterator (1.2x slower on
    // d=2 18SV9) variants regressed performance versus the raw loop.
  }


  auto build_remaining_amplicons_list(uint64_t const swarmed,
                                      uint64_t const seed_abundance,
                                      struct Parameters const & parameters,
                                      Data const & data,
                                      std::vector<struct ampliconinfo_s> const & amps_v,
                                      Cluster_workspace & workspace) -> uint64_t {
    auto const amplicons = amps_v.size();
    uint64_t listlen {0};
    for (auto i = swarmed; i < amplicons; ++i) {
      auto const ampliconid = amps_v[i].ampliconid;
      if (parameters.opt_no_cluster_breaking or
          (data.abundance(ampliconid) <= seed_abundance)) {
        workspace.qgramamps_v[listlen] = ampliconid;
        ++listlen;
      }
    }
    return listlen;
  }


  auto build_subseed_candidate_list(uint64_t const swarmed,
                                    uint64_t const amplicons,
                                    struct ampliconinfo_s const & subseed,
                                    struct Parameters const & parameters,
                                    Data const & data,
                                    std::vector<struct ampliconinfo_s> const & amps_v,
                                    Cluster_workspace & workspace) -> uint64_t {
    auto const subseed_abundance = data.abundance(subseed.ampliconid);
    uint64_t subseedlistlen {0};
    for (auto i = swarmed; i < amplicons; ++i) {
      uint64_t const targetampliconid = amps_v[i].ampliconid;
      if ((amps_v[i].diffestimate <=
           subseed.radius + parameters.opt_differences) and
          (parameters.opt_no_cluster_breaking or
           (data.abundance(targetampliconid)
            <= subseed_abundance))) {
        workspace.qgramamps_v[subseedlistlen] = targetampliconid;
        workspace.qgramindices_v[subseedlistlen] = i;
        ++subseedlistlen;
      }
    }
    return subseedlistlen;
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


  auto finalize_algo_run(uint64_t const amplicons,
                         unsigned int const swarmid,
                         uint64_t const largestswarm,
                         uint64_t const maxgenerations,
                         struct Parameters const & parameters,
                         Data const & data,
                         std::vector<struct ampliconinfo_s> const & amps_v) -> void {
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

    static_cast<void>(std::fputc('\n', parameters.logfile));

    std::fprintf(parameters.logfile, "Number of swarms:  %u\n", swarmid);

    std::fprintf(parameters.logfile, "Largest swarm:     %" PRIu64 "\n", largestswarm);

    std::fprintf(parameters.logfile, "Max generations:   %" PRIu64 "\n", maxgenerations);
  }


  auto start_new_cluster(Pool_cursor & cursor,
                         unsigned int const swarmid,
                         std::vector<struct ampliconinfo_s> & amps_v,
                         Cluster_state & state,
                         Cluster_workspace & workspace,
                         Data const & data) -> uint64_t {
    /* process each initial seed */
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


  auto seed_first_generation(struct Parameters const & parameters,
                             Data const & data,
                             QgramDiffer & qgram_differ,
                             Scanner & scanner,
                             Bit_mode const bits,
                             Pool_cursor & cursor,
                             unsigned int const swarmid,
                             uint64_t const seedindex,
                             std::vector<struct ampliconinfo_s> & amps_v,
                             Cluster_workspace & workspace,
                             Cluster_state & state) -> void {
    /* find diff estimates between seed and each amplicon in pool */
    uint64_t const seedampliconid = amps_v[seedindex].ampliconid;
    auto const seed_abundance = data.abundance(seedampliconid);

    uint64_t const listlen = build_remaining_amplicons_list(cursor.swarmed, seed_abundance,
                                                            parameters, data,
                                                            amps_v, workspace);

    // the collected candidates, not the pool-sized scratch behind them
    qgram_differ.fast(seedampliconid,
                      View<uint64_t>{workspace.qgramamps_v.data(), listlen},
                      Span<uint64_t>{workspace.qgramdiffs_v.data(), listlen});

    uint64_t targetcount = 0;
    for (auto i = 0ULL; i < listlen; ++i) {
      auto const poolampliconid = workspace.qgramamps_v[i];
      auto const diff = workspace.qgramdiffs_v[i];
      assert(diff <= std::numeric_limits<unsigned int>::max());
      amps_v[cursor.swarmed + i].diffestimate = static_cast<unsigned int>(diff);
      if (diff <= parameters.opt_differences) {
        workspace.targetindices[targetcount] = cursor.swarmed + i;
        workspace.targetampliconids[targetcount] = poolampliconid;
        ++targetcount;
      }
    }

    if (targetcount == 0) { return; }

    // the candidate window inside the workspace: the first targetcount
    // entries of the target list and of its three result columns
    auto const targets = View<uint64_t>{workspace.targetampliconids.data(), targetcount};
    auto const scores = Span<uint64_t>{workspace.scores_v.data(), targetcount};
    auto const diffs = Span<uint64_t>{workspace.diffs_v.data(), targetcount};
    auto const alignlengths = Span<uint64_t>{workspace.alignlengths.data(), targetcount};
    scanner.run(seedampliconid, targets, scores, diffs, alignlengths, bits);

    for (auto target_id = 0ULL; target_id < targetcount; ++target_id) {
      auto const diff = workspace.diffs_v[target_id];

      if (diff > parameters.opt_differences) { continue; }
      auto const target = workspace.targetindices[target_id];

      /* move the 'target' to the position ('swarmed')
         of the first unswarmed amplicon in the pool */
      move_target_to_first_unswarmed_position(cursor.swarmed, target, amps_v);

      include_amplicon_in_cluster(cursor.swarmed, diff, swarmid,
                                  amps_v[seedindex], amps_v, workspace.hits,
                                  state, parameters, data);
      ++cursor.swarmed;
    }
  }


  auto grow_cluster_from_subseeds(struct Parameters const & parameters,
                                  Data const & data,
                                  QgramDiffer & qgram_differ,
                                  Scanner & scanner,
                                  Bit_mode const bits,
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
                                                               subseed, parameters,
                                                               data, amps_v, workspace);

      // as above: the collected candidates, not the whole scratch buffer
      qgram_differ.fast(subseed.ampliconid,
                        View<uint64_t>{workspace.qgramamps_v.data(), subseedlistlen},
                        Span<uint64_t>{workspace.qgramdiffs_v.data(), subseedlistlen});

      for (auto i = 0ULL; i < subseedlistlen; ++i) {
        if (workspace.qgramdiffs_v[i] <= parameters.opt_differences) {
          workspace.targetindices[targetcount] = workspace.qgramindices_v[i];
          workspace.targetampliconids[targetcount] = workspace.qgramamps_v[i];
          ++targetcount;
        }
      }

      if (targetcount == 0) { continue; }

      // the candidate window inside the workspace, as above
      auto const targets = View<uint64_t>{workspace.targetampliconids.data(), targetcount};
      auto const scores = Span<uint64_t>{workspace.scores_v.data(), targetcount};
      auto const diffs = Span<uint64_t>{workspace.diffs_v.data(), targetcount};
      auto const alignlengths = Span<uint64_t>{workspace.alignlengths.data(), targetcount};
      scanner.run(subseed.ampliconid, targets, scores, diffs, alignlengths, bits);

      for (auto target_id = 0ULL; target_id < targetcount; ++target_id) {
        auto const diff = workspace.diffs_v[target_id];

        if (diff > parameters.opt_differences) { continue; }
        auto const target = workspace.targetindices[target_id];

        /* find correct position in list */

        auto const pos = find_correct_position_in_list(cursor.swarmed, target, cursor.seeded,
                                                       subseed, amps_v);

        move_target_to_first_unswarmed_position(pos, target, amps_v);

        include_amplicon_in_cluster(pos, diff, swarmid, subseed,
                                    amps_v, workspace.hits, state,
                                    parameters, data);
        ++cursor.swarmed;
      }
    }
  }


} // namespace


auto algo_run(struct Parameters const & parameters,
              Data const & data) -> void {
  // RAII: allocates per-thread scratch buffers and starts the worker
  // pool; threads are joined when the Scanner is destroyed at end of
  // algo_run scope.
  Scanner scanner(parameters, data);

  uint64_t largestswarm {0};
  uint64_t maxgenerations {0};

  uint64_t const amplicons = data.sequence_count();
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
    // NwAligner's direction matrix is O(longestamplicon^2); guard it like
    // the search buffers so a pathological sequence length fails with a
    // clear message instead of aborting inside operator new.
    require_ram(longestamplicon * longestamplicon, 1,
               "the UCLUST alignment matrix");
    aligner = utils::make_unique<NwAligner>(
        longestamplicon,
        parameters.penalty_mismatch,
        static_cast<unsigned long int>(parameters.penalty_gapopen),
        static_cast<unsigned long int>(parameters.penalty_gapextend));
  }

  set_amplicon_ids(amps_v);
  auto const bits = set_bit_mode(parameters);

  Pool_cursor cursor;

  auto swarmid = 0U;

  Progress progress("Clustering:       ", amplicons, parameters);
  while (cursor.seeded < amplicons) {

      ++swarmid;
      Cluster_state state;
      uint64_t const seedindex = cursor.seeded;

      uint64_t const seedampliconid = start_new_cluster(cursor, swarmid, amps_v,
                                                        state, workspace, data);

      seed_first_generation(parameters, data, qgram_differ,
                            scanner, bits,
                            cursor, swarmid, seedindex,
                            amps_v, workspace, state);

      grow_cluster_from_subseeds(parameters, data, qgram_differ,
                                 scanner, bits,
                                 cursor, swarmid, amplicons,
                                 amps_v, workspace, state);

      largestswarm = std::max(state.swarmsize, largestswarm);
      maxgenerations = std::max(state.maxgen, maxgenerations);

      write_cluster_outputs(swarmid, seedampliconid, state, workspace.hits,
                            aligner.get(), parameters, data);

      progress.update(cursor.seeded);
  }
  progress.done();

  finalize_algo_run(amplicons, swarmid, largestswarm, maxgenerations,
                    parameters, data, amps_v);
}

