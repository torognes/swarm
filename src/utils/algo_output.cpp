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

#include "algo_output.hpp"
#include "algo_internal.hpp"
#include "../swarm.hpp"
#include "../db.hpp"
#include "nw_aligner.hpp"
#include "print_view.hpp"  // fprint
#include "progress.hpp"
#include "view.hpp"  // View, make_view
#include <algorithm>  // std::sort
#include <cassert>
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fprintf()
#include <vector>


namespace {

  auto collect_seeds(struct Parameters const & parameters,
                     Data const & data,
                     View<struct ampliconinfo_s> const amps) -> std::vector<struct swarminfo_t> {
    auto const amplicons = amps.size();
    Progress progress("Collecting seeds:    ", amplicons, parameters);
    std::vector<struct swarminfo_t> seeds(amplicons);
    auto swarmcount = 0UL;
    uint64_t mass = 0;
    auto previous_id = amps.front().swarmid;
    auto seed = amps.front().ampliconid;
    mass += data.abundance(seed);
    // the first amplicon is accounted for above, so the loop covers the
    // rest. progress.increment() replaces update(i): it counts from one on
    // its first call, which is where the index loop started.
    for (auto const & amplicon : amps.drop(1)) {
        auto const current_id = amplicon.swarmid;
        if (current_id != previous_id) {
            seeds[swarmcount].seed = seed;  // update previous
            seeds[swarmcount].mass = mass;
            ++swarmcount;
            mass = 0;
            seed = amplicon.ampliconid;
          }
        mass += data.abundance(amplicon.ampliconid);
        previous_id = current_id;
        progress.increment();
      }
    seeds[swarmcount].seed = seed;
    seeds[swarmcount].mass = mass;
    ++swarmcount;

    // free some memory
    assert(swarmcount <= seeds.size());  // resize() would lengthen, not shrink
    seeds.resize(swarmcount);
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
    auto * const seeds_file = parameters.seeds_file.get();
    // one scratch buffer for the whole file, as with NwAligner below
    Sequence_printer const sequence_printer {data.longest_sequence()};

    for (auto const& seed: seeds) {
      auto const swarm_mass = seed.mass;
      auto const swarm_seed = seed.seed;

      fprint(seeds_file, '>');
      fprint_id_with_new_abundance(seeds_file, data.info(swarm_seed), swarm_mass,
                                   parameters.opt_usearch_abundance);
      fprint(seeds_file, '\n');
      sequence_printer.print(seeds_file, data.sequence_view(swarm_seed));
      progress.increment();
    }
    progress.done();
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
    auto * const stats_file = parameters.statsfile.get();

    fprint_integer(stats_file, swarmsize);
    fprint(stats_file, '\t');
    fprint_integer(stats_file, amplicons_copies);
    fprint(stats_file, '\t');
    fprint_id_noabundance(stats_file, data.info(seedampliconid), parameters.opt_usearch_abundance);
    fprint(stats_file, '\t');
    fprint_integer(stats_file, abundance);
    fprint(stats_file, '\t');
    fprint_integer(stats_file, singletons);
    fprint(stats_file, '\t');
    fprint_integer(stats_file, maxgen);
    fprint(stats_file, '\t');
    fprint_integer(stats_file, maxradius);
    fprint(stats_file, '\n');
  }


  auto write_uclust_cluster(unsigned int const swarmid,
                            uint64_t const swarmsize,
                            uint64_t const seedampliconid,
                            View<uint64_t> const hits,
                            NwAligner & aligner,
                            struct Parameters const & parameters,
                            Data const & data) -> void {
    auto const seed_seq = data.sequence_view(seedampliconid);
    auto * const uclust_file = parameters.uclustfile.get();
    // uclust numbers its clusters from zero, swarmid from one
    auto const cluster_no = swarmid - 1;

    fprint(uclust_file, "C\t");
    fprint_integer(uclust_file, cluster_no);
    fprint(uclust_file, '\t');
    fprint_integer(uclust_file, swarmsize);
    fprint(uclust_file, "\t*\t*\t*\t*\t*\t");
    fprint_id(uclust_file, data.info(seedampliconid), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    fprint(uclust_file, "\t*\n");

    fprint(uclust_file, "S\t");
    fprint_integer(uclust_file, cluster_no);
    fprint(uclust_file, '\t');
    fprint_integer(uclust_file, seed_seq.length);
    fprint(uclust_file, "\t*\t*\t*\t*\t*\t");
    fprint_id(uclust_file, data.info(seedampliconid), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    fprint(uclust_file, "\t*\n");

    // the cluster's members except its seed, which the S line above
    // already reported
    for (auto const hit : hits.drop(1)) {
      auto const hit_seq = data.sequence_view(hit);

      auto const result = aligner.align(hit_seq, seed_seq);

      fprint(uclust_file, "H\t");
      fprint_integer(uclust_file, cluster_no);
      fprint(uclust_file, '\t');
      fprint_integer(uclust_file, hit_seq.length);
      fprint(uclust_file, '\t');
      // the one field fprintf still formats here, and the one thing it does
      // that nothing simpler does: a double, rounded to one decimal. Doing
      // it by hand would mean reproducing glibc's rounding, for a field
      // written once per cluster member right after a full alignment.
      static_cast<void>(std::fprintf(uclust_file, "%.1f", result.percent_id));
      fprint(uclust_file, "\t+\t0\t0\t");
      // a hit is never identical to its centroid: dereplication is
      // mandatory at d >= 1 (db.cpp, "some fasta entries have identical
      // sequences"), and it compares packed encodings, so U/T and case
      // variants count as duplicates too. USEARCH's '=' CIGAR, which marks
      // an identical hit, is therefore emitted only by dereplicate.cpp.
      assert(result.differences > 0);
      fprint(uclust_file, result.cigar_string);
      fprint(uclust_file, '\t');

      fprint_id(uclust_file, data.info(hit), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      fprint(uclust_file, '\t');
      fprint_id(uclust_file, data.info(seedampliconid), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      fprint(uclust_file, '\n');
    }
  }


  // Separators in a one-line-per-swarm listing: 'within' joins amplicons
  // of the same swarm, 'between' starts the next swarm.
  struct Swarm_separators {
    char within;
    char between;
  };


  auto write_swarm_listing(Swarm_separators const separators,
                           struct Parameters const & parameters,
                           Data const & data,
                           View<struct ampliconinfo_s> const amps) -> void {
    fprint_id(parameters.outfile.get(), data.info(amps.front().ampliconid),
              parameters.opt_usearch_abundance, parameters.opt_append_abundance);
    auto previous_id = amps.front().swarmid;

    // the first amplicon is printed above, so the loop covers the rest
    for (auto const & amplicon : amps.drop(1)) {
        auto const current_id = amplicon.swarmid;
        fprint(parameters.outfile.get(),
               current_id == previous_id ? separators.within : separators.between);
        fprint_id(parameters.outfile.get(), data.info(amplicon.ampliconid),
                       parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        previous_id = current_id;
      }
    fprint(parameters.outfile.get(), '\n');
  }

} // namespace


auto write_swarms_default_format(struct Parameters const & parameters,
                                 Data const & data,
                                 View<struct ampliconinfo_s> const amps) -> void {
  /* native swarm output */
  static constexpr Swarm_separators separators {' ' /* usually a space */, '\n'};
  write_swarm_listing(separators, parameters, data, amps);
}


auto write_swarms_mothur_format(unsigned int const swarmid,
                                struct Parameters const & parameters,
                                Data const & data,
                                View<struct ampliconinfo_s> const amps) -> void {
  /* mothur list file output */
  static constexpr Swarm_separators separators {',', '\t'};
  fprint(parameters.outfile.get(), "swarm_");
  fprint_integer(parameters.outfile.get(), parameters.opt_differences);
  fprint(parameters.outfile.get(), '\t');
  fprint_integer(parameters.outfile.get(), swarmid);
  fprint(parameters.outfile.get(), '\t');
  write_swarm_listing(separators, parameters, data, amps);
}


auto write_internal_structure_line(uint64_t const parent_id,
                                   uint64_t const child_id,
                                   uint64_t const diff,
                                   unsigned int const swarmid,
                                   unsigned int const generation,
                                   struct Parameters const & parameters,
                                   Data const & data) -> void {
  // hoisted: the handle is named eleven times below, and .get() on a
  // unique_ptr that is not reassigned here yields the same pointer each time
  auto * const structure_file = parameters.internal_structure_file.get();

  fprint_id_noabundance(structure_file, data.info(parent_id), parameters.opt_usearch_abundance);
  fprint(structure_file, '\t');
  fprint_id_noabundance(structure_file, data.info(child_id), parameters.opt_usearch_abundance);
  fprint(structure_file, '\t');
  fprint_integer(structure_file, diff);
  fprint(structure_file, '\t');
  fprint_integer(structure_file, swarmid);
  fprint(structure_file, '\t');
  fprint_integer(structure_file, generation);
  fprint(structure_file, '\n');
}


auto write_representative_sequences(struct Parameters const & parameters,
                                    Data const & data,
                                    View<struct ampliconinfo_s> const amps) -> void {
  auto seeds = collect_seeds(parameters, data, amps);
  sort_seeds(parameters, data, seeds);
  write_seeds(parameters, data, seeds);
}


auto write_cluster_outputs(unsigned int const swarmid,
                           uint64_t const seedampliconid,
                           Cluster_state const & state,
                           View<uint64_t> const hits,
                           NwAligner * const aligner,
                           struct Parameters const & parameters,
                           Data const & data) -> void {
  if (parameters.uclustfile.get() != nullptr) {
    write_uclust_cluster(swarmid, state.swarmsize, seedampliconid, hits,
                         *aligner, parameters, data);
  }

  if (parameters.statsfile.get() != nullptr) {
    write_stats_line(state.swarmsize, state.amplicons_copies, state.singletons,
                     state.maxgen, state.maxradius, seedampliconid, parameters, data);
  }
}
