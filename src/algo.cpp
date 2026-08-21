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
#include "utils/search_data.hpp"  // score_ceiling_8
#include "utils/memory_budget.hpp"
#include "utils/make_unique.hpp"
#include "utils/print_view.hpp"  // fprint, fprint_integer
#include "utils/progress.hpp"
#include "utils/span.hpp"  // Span, make_span
#include "utils/view.hpp"  // View, make_view
#include <algorithm>  // std::min(), std::for_each, std::is_sorted, std::partition_point
#include <cassert>
#include <cstdint>  // int64_t, uint64_t
#include <iterator>  // std::next, std::prev, std::distance
#include <memory>  // unique pointer
#include <string>
#include <vector>

#ifndef NDEBUG
#include <limits>  // std::numeric_limits, in the assertions only
#endif


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
    // A lane that saturates sticks at exactly score_ceiling_8, so a genuine
    // score of score_ceiling_8 would be indistinguishable from an overflow
    // (see save_score_8). Cap the worst-case score at score_ceiling_8 - 1 so
    // that a score of score_ceiling_8 can only ever mean "saturated".
    static constexpr auto max_reliable_score = score_ceiling_8 - 1;

    // The guard below makes a release build on aarch64 disagree with a debug
    // or coverage build on the same platform, and that is deliberate rather
    // than a leftover. A release build takes the 16-bit kernel for every d,
    // which is the one aarch64 users run; a debug or coverage build falls
    // through to the saturation rule below and takes the 8-bit kernel at low
    // d, which is how that kernel stays exercised on this platform. Between
    // the two builds the suite covers both.
    //
    // So do not "fix" the divergence by dropping the guard: that would leave
    // the 8-bit kernel untested on aarch64. Anyone reasoning about which
    // kernel a given aarch64 run used has to look at how it was built.
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
      = static_cast<uint64_t>(std::min(max_reliable_score / parameters.penalty_mismatch,
                                       max_reliable_score / (parameters.penalty_gapopen +
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


#ifndef NDEBUG
  // Whether the whole amplicon pool is in decreasing abundance order, which
  // sort_index_if_need_be() (db.cpp) establishes and the clustering loop
  // preserves. This is the premise of first_qualifying_amplicon() below and
  // of the seed's candidate list in seed_first_generation(), and it is
  // checked once for the whole run rather than at either use: both of those
  // run once per subseed or per cluster, so an O(pool) check there made a
  // debug build quadratic on top of the walks it was guarding.
  auto pool_is_abundance_sorted(Data const & data,
                                std::vector<struct ampliconinfo_s> const & amps_v) -> bool {
    return std::is_sorted(amps_v.cbegin(), amps_v.cend(),
                          [&data](struct ampliconinfo_s const & lhs,
                                  struct ampliconinfo_s const & rhs) -> bool {
                            return data.abundance(lhs.ampliconid)
                                 > data.abundance(rhs.ampliconid);
                          });
  }
#endif


  // The position in the pool of the first target that cluster breaking
  // admits for this subseed, i.e. the first entry no more abundant than the
  // subseed itself, or the end of the pool if there is none. The pool is in
  // decreasing abundance order (see pool_is_abundance_sorted), so those
  // targets are a suffix of it and binary search finds where it opens.
  auto first_qualifying_amplicon(uint64_t const swarmed,
                                 struct ampliconinfo_s const & subseed,
                                 Data const & data,
                                 std::vector<struct ampliconinfo_s> const & amps_v) -> uint64_t {
    auto const subseed_abundance = data.abundance(subseed.ampliconid);
    auto const pool_start = std::next(amps_v.cbegin(),
                                      static_cast<std::ptrdiff_t>(swarmed));
    auto const suffix =
      std::partition_point(pool_start, amps_v.cend(),
                           [&data, subseed_abundance](struct ampliconinfo_s const & amplicon) -> bool {
                             return data.abundance(amplicon.ampliconid) > subseed_abundance;
                           });

    // The two entries either side of the boundary, which is as much as can be
    // checked here without walking the range std::partition_point just
    // bisected. They do not hold for an unpartitioned range except by
    // coincidence, so between them and the once-per-run ordering check a
    // broken premise has to get past both to go unnoticed.
    assert((suffix == amps_v.cend()) or
           (data.abundance(suffix->ampliconid) <= subseed_abundance));
    assert((suffix == pool_start) or
           (data.abundance(std::prev(suffix)->ampliconid) > subseed_abundance));

    return static_cast<uint64_t>(std::distance(amps_v.cbegin(), suffix));
  }


  auto build_subseed_candidate_list(uint64_t const swarmed,
                                    uint64_t const amplicons,
                                    struct ampliconinfo_s const & subseed,
                                    struct Parameters const & parameters,
                                    Data const & data,
                                    std::vector<struct ampliconinfo_s> const & amps_v,
                                    Cluster_workspace & workspace) -> uint64_t {
    // Where the walk can start. Testing the abundance of every entry was
    // 2.1e9 random loads over a -d 2 run on 219k reads, and the entries in
    // front of the suffix -- 28 % of the walk -- were stepped over only to
    // be rejected. Both go away: inside the suffix the test is redundant,
    // and opt_no_cluster_breaking makes it vacuous everywhere.
    auto const first = parameters.opt_no_cluster_breaking
      ? swarmed
      : first_qualifying_amplicon(swarmed, subseed, data, amps_v);

    uint64_t subseedlistlen {0};
    for (auto i = first; i < amplicons; ++i) {
      if (amps_v[i].diffestimate <=
          subseed.radius + parameters.opt_differences) {
        workspace.qgramamps_v[subseedlistlen] = amps_v[i].ampliconid;
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


  auto finalize_algo_run(unsigned int const swarmid,
                         uint64_t const largestswarm,
                         uint64_t const maxgenerations,
                         struct Parameters const & parameters,
                         Data const & data,
                         View<struct ampliconinfo_s> const amps) -> void {
    /* output swarms */
    if (not amps.empty()) {
      if (parameters.opt_mothur) {
        write_swarms_mothur_format(swarmid, parameters, data, amps);
      }
      else {
        write_swarms_default_format(parameters, data, amps);
      }
    }


    /* dump seeds in fasta format with sum of abundances */
    if ((not parameters.opt_seeds.empty()) and (not amps.empty())) {
      write_representative_sequences(parameters, data, amps);
    }

    fprint(parameters.logfile, '\n');

    fprint(parameters.logfile, "Number of swarms:  ");
    fprint_integer(parameters.logfile, swarmid);
    fprint(parameters.logfile, '\n');

    fprint(parameters.logfile, "Largest swarm:     ");
    fprint_integer(parameters.logfile, largestswarm);
    fprint(parameters.logfile, '\n');

    fprint(parameters.logfile, "Max generations:   ");
    fprint_integer(parameters.logfile, maxgenerations);
    fprint(parameters.logfile, '\n');
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

    // The seed is at least as abundant as every entry of the pool it is
    // about to be compared against, which is why the candidate list below is
    // the whole of that pool: the abundance test that used to select
    // candidates here rejected nothing -- 0 of 4.7e9 entries over a d = 2 run
    // on 219k reads -- and opt_no_cluster_breaking never entered into it, it
    // only widened a filter that was already total.
    //
    // Testing the pool front is enough because the pool is in decreasing
    // abundance order, which pool_is_abundance_sorted() asserts once for the
    // whole run. Three facts keep it that way from there on:
    // sort_index_if_need_be() (db.cpp) always leaves the amplicons in
    // decreasing abundance order, move_target_to_first_unswarmed_position()
    // lifts one entry out and shifts the rest, and dropping one element from
    // an ordered sequence leaves the rest ordered.
    assert((cursor.swarmed == amps_v.size()) or
           (data.abundance(amps_v[cursor.swarmed].ampliconid)
            <= data.abundance(seedampliconid)));

    // The seed's candidates are the whole unswarmed pool, in pool order, so
    // there is nothing to collect: hand the workers the pool itself. The pass
    // this replaces copied four bytes out of each twenty-byte record into a
    // contiguous list -- 4 711 615 447 iterations over a -d 2 run on 219k
    // reads -- and did it on this thread, while the workers were about to
    // read those records' q-gram vectors anyway.
    //
    // It also puts the pairing beyond doubt. The loop below writes
    // amps_v[cursor.swarmed + i].diffestimate for candidate i, which is that
    // candidate's own entry because the view handed to the scan *is* the pool
    // from cursor.swarmed on. That used to rest on the collecting pass
    // keeping every entry, asserted here, and now holds by construction --
    // had a candidate ever been left out, every diffestimate past it would
    // have landed on the wrong amplicon and the subseed passes, which prune
    // on diffestimate, would have dropped true neighbours.
    auto const listlen = amps_v.size() - cursor.swarmed;
    qgram_differ.fast_over_pool(seedampliconid,
                                make_view(amps_v).subview(cursor.swarmed, listlen),
                                workspace.qgram_diffs(listlen));

    uint64_t targetcount = 0;
    for (auto i = 0ULL; i < listlen; ++i) {
      auto const diff = workspace.qgramdiffs_v[i];
      assert(diff <= std::numeric_limits<unsigned int>::max());
      // the id comes from the record this iteration is writing to, rather
      // than from a second array the collecting pass used to fill
      auto & pool_entry = amps_v[cursor.swarmed + i];
      pool_entry.diffestimate = static_cast<unsigned int>(diff);
      if (diff <= parameters.opt_differences) {
        workspace.targetindices[targetcount] = cursor.swarmed + i;
        workspace.targetampliconids[targetcount] = pool_entry.ampliconid;
        ++targetcount;
      }
    }

    if (targetcount == 0) { return; }

    // the candidate window inside the workspace: the first targetcount
    // entries of the target list and of its three result columns
    scanner.run(seedampliconid,
                workspace.targets(targetcount),
                workspace.scores(targetcount),
                workspace.diffs(targetcount),
                bits);

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
                        workspace.qgram_candidates(subseedlistlen),
                        workspace.qgram_diffs(subseedlistlen));

      for (auto i = 0ULL; i < subseedlistlen; ++i) {
        if (workspace.qgramdiffs_v[i] <= parameters.opt_differences) {
          workspace.targetindices[targetcount] = workspace.qgramindices_v[i];
          workspace.targetampliconids[targetcount] = workspace.qgramamps_v[i];
          ++targetcount;
        }
      }

      if (targetcount == 0) { continue; }

      // the candidate window inside the workspace, as above
      scanner.run(subseed.ampliconid,
                  workspace.targets(targetcount),
                  workspace.scores(targetcount),
                  workspace.diffs(targetcount),
                  bits);

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
  assert(pool_is_abundance_sorted(data, amps_v));
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

      // the cluster's members, seed first: the filled part of the hit
      // buffer, not the pool-sized buffer behind it
      write_cluster_outputs(swarmid, seedampliconid, state,
                            make_view(workspace.hits).first(state.hitcount),
                            aligner.get(), parameters, data);

      progress.update(cursor.seeded);
  }
  progress.done();

  finalize_algo_run(swarmid, largestswarm, maxgenerations,
                    parameters, data, make_view(amps_v));
}

