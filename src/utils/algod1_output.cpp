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

#include "algod1_output.hpp"
#include "../swarm.hpp"
#include "../db.hpp"
#include "algod1_internal.hpp"
#include "nw_aligner.hpp"
#include "print_view.hpp"  // fprint
#include "progress.hpp"
#include "span.hpp"
#include <algorithm>  // std::sort()
#include <cassert>  // assert()
#include <cstdint>  // uint64_t
#include <cstdio>  // fprintf()
#include <numeric>  // std::iota
#include <vector>


auto write_network_file(const uint64_t number_of_networks,
                        struct Parameters const & parameters,
                        Data const & data,
                        std::vector<struct ampinfo_s> const & ampinfo_v,
                        std::vector<unsigned int> & network_v) -> void {
  // a network is a cluster with at least two sequences (no singletons)
  Progress progress("Dumping network:  ", number_of_networks, parameters);

  assert(ampinfo_v.size() == data.sequence_count());
  auto counter = 0ULL;
  for (auto const& amplicon: ampinfo_v) {
    auto const neighbours = neighbours_of(network_v, amplicon);

    // amplicon indexes are already sorted by decreasing abundance
    // then by header in db.cpp, so a natural ascending sort here
    // emits neighbours in that ranking order. Earlier dereplication
    // guarantees indexes are distinct.
    std::sort(neighbours.begin(), neighbours.end());

    for (auto const neighbour : neighbours)
      {
        fprint_id(parameters.network_file.get(), data.info(counter), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        fprint(parameters.network_file.get(), '\t');
        fprint_id(parameters.network_file.get(), data.info(neighbour), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        fprint(parameters.network_file.get(), '\n');
        progress.increment();
      }
    ++counter;
  }
  progress.done();
}


namespace {

  // swarminfo_v is indexed by ampinfo_s::swarmid, an unsigned int that
  // reserves no_swarm, so a cluster count always fits in an unsigned int.
  // The writers below take their bound from here rather than from
  // swarminfo_v.size(): that keeps each loop counter the same type as the
  // ids it hands to cluster_members(), instead of walking an unsigned int
  // towards a 64-bit size_type it could never reach.
  auto cluster_count(std::vector<struct swarminfo_s> const & swarminfo_v) noexcept -> unsigned int {
    assert(swarminfo_v.size() < no_swarm);
    return static_cast<unsigned int>(swarminfo_v.size());
  }


  auto write_swarms_default_format(struct Parameters const & parameters,
                                   Data const & data,
                                   std::vector<struct ampinfo_s> const & ampinfo_v,
                                   std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    static constexpr char sepchar {' '};
    auto const clusters = cluster_count(swarminfo_v);
    Progress progress("Writing swarms:   ", clusters, parameters);

    for (auto swarmid = 0U; swarmid < clusters; ++swarmid) {
      if (swarminfo_v[swarmid].attached) {
        continue;
      }

      const auto seed = swarminfo_v[swarmid].seed;
      for (auto const amp_id : cluster_members(ampinfo_v, seed)) {
        if (amp_id != seed) {
          fprint(parameters.outfile.get(), sepchar);
        }
        fprint_id(parameters.outfile.get(), data.info(amp_id),
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      }
      fprint(parameters.outfile.get(), '\n');
      progress.update(swarmid + 1);
    }

    progress.done();
  }


  auto write_swarms_mothur_format(struct Parameters const & parameters,
                                  Data const & data,
                                  std::vector<struct ampinfo_s> const & ampinfo_v,
                                  std::vector<struct swarminfo_s> const & swarminfo_v,
                                  Overall_stats const & overall_stats) -> void {
    auto const clusters = cluster_count(swarminfo_v);
    Progress progress("Writing swarms:   ", clusters, parameters);

    fprint(parameters.outfile.get(), "swarm_");
    fprint_integer(parameters.outfile.get(), parameters.opt_differences);
    fprint(parameters.outfile.get(), '\t');
    fprint_integer(parameters.outfile.get(), overall_stats.swarmcount_adjusted);

    for (auto swarmid = 0U; swarmid < clusters; ++swarmid) {
      assert(not swarminfo_v[swarmid].attached);
      if (swarminfo_v[swarmid].attached) {
        continue;
      }

      const auto seed = swarminfo_v[swarmid].seed;
      for (auto const amp_id : cluster_members(ampinfo_v, seed)) {
        if (amp_id == seed) {
          fprint(parameters.outfile.get(), '\t');
        }
        else {
          fprint(parameters.outfile.get(), ',');
        }
        fprint_id(parameters.outfile.get(), data.info(amp_id),
                  parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      }
      progress.update(swarmid + 1);
    }

    fprint(parameters.outfile.get(), '\n');

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
    auto * const uclust_file = parameters.uclustfile.get();

    for (auto const & swarm_info : swarminfo_v) {
      if (swarm_info.attached) {
        continue;
      }

      const auto seed = swarm_info.seed;

      fprint(uclust_file, "C\t");
      fprint_integer(uclust_file, cluster_no);
      fprint(uclust_file, '\t');
      fprint_integer(uclust_file, swarm_info.size);
      fprint(uclust_file, "\t*\t*\t*\t*\t*\t");
      fprint_id(uclust_file, data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      fprint(uclust_file, "\t*\n");

      fprint(uclust_file, "S\t");
      fprint_integer(uclust_file, cluster_no);
      fprint(uclust_file, '\t');
      fprint_integer(uclust_file, data.sequence_view(seed).length);
      fprint(uclust_file, "\t*\t*\t*\t*\t*\t");
      fprint_id(uclust_file, data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      fprint(uclust_file, "\t*\n");

      auto const seed_seq = data.sequence_view(seed);
      for (auto const amp_id : cluster_members_after_seed(ampinfo_v, seed))
        {
          auto const amp_seq = data.sequence_view(amp_id);

          auto const result = aligner.align(amp_seq, seed_seq);

          fprint(uclust_file, "H\t");
          fprint_integer(uclust_file, cluster_no);
          fprint(uclust_file, '\t');
          fprint_integer(uclust_file, amp_seq.length);
          fprint(uclust_file, '\t');
          // the one field fprintf still formats: a double to one decimal.
          // See the identical note in algo_output.cpp.
          static_cast<void>(std::fprintf(uclust_file, "%.1f", result.percent_id));
          fprint(uclust_file, "\t+\t0\t0\t");
          if (result.differences > 0) {
            fprint(uclust_file, result.cigar_string);
          }
          else {
            fprint(uclust_file, '=');
          }
          fprint(uclust_file, '\t');

          fprint_id(uclust_file, data.info(amp_id), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          fprint(uclust_file, '\t');
          fprint_id(uclust_file, data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          fprint(uclust_file, '\n');
        }

      ++cluster_no;
      progress.increment();
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

    auto * const seeds_file = parameters.seeds_file.get();
    // one scratch buffer for the whole file
    Sequence_printer const sequence_printer {data.longest_sequence()};

    for (const auto index : sorter) {
      const auto & a_swarm = swarminfo_v[index];
      if (a_swarm.attached) {
        continue;
      }
      const auto seed = a_swarm.seed;
      const auto mass = a_swarm.mass;
      fprint(seeds_file, '>');
      fprint_id_with_new_abundance(seeds_file, data.info(seed), mass,
                                   parameters.opt_usearch_abundance);
      fprint(seeds_file, '\n');
      sequence_printer.print(seeds_file, data.sequence_view(seed));
      progress.increment();
    }

    progress.done();
  }


  auto write_structure_file(struct Parameters const & parameters,
                            Data const & data,
                            std::vector<struct ampinfo_s> const & ampinfo_v,
                            std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    auto cluster_no = 0U;

    auto const clusters = cluster_count(swarminfo_v);
    Progress progress("Writing structure:", clusters, parameters);

    // Column 3 of the internal-structure file is the number of differences
    // between the two amplicons (man swarm, --internal-structure). It is a
    // constant per link kind rather than a computed value: a fastidious
    // graft joins amplicons two differences apart, an ordinary link one.
    // Named because "1" and "2" in adjacent printf arguments are exactly
    // the pair a reader cannot tell apart.
    static constexpr unsigned int graft_differences {2};
    static constexpr unsigned int link_differences {1};

    auto * const structure_file = parameters.internal_structure_file.get();

    for (auto swarmid = 0U; swarmid < clusters; ++swarmid)
      {
        if (swarminfo_v[swarmid].attached) {
          continue;
        }
        const auto seed = swarminfo_v[swarmid].seed;

        for (auto const amp_id : cluster_members_after_seed(ampinfo_v, seed))
          {
            const auto graft_parent = ampinfo_v[amp_id].graft_cand;
            if (graft_parent != no_swarm)
              {
                fprint_id_noabundance(structure_file, data.info(graft_parent),
                                           parameters.opt_usearch_abundance);
                fprint(structure_file, '\t');
                fprint_id_noabundance(structure_file, data.info(amp_id), parameters.opt_usearch_abundance);
                fprint(structure_file, '\t');
                fprint_integer(structure_file, graft_differences);
                fprint(structure_file, '\t');
                fprint_integer(structure_file, cluster_no + 1);
                fprint(structure_file, '\t');
                fprint_integer(structure_file, ampinfo_v[graft_parent].generation + 1);
                fprint(structure_file, '\n');
              }

            const auto parent = ampinfo_v[amp_id].parent;
            if (parent != no_swarm)
              {
                fprint_id_noabundance(structure_file, data.info(parent), parameters.opt_usearch_abundance);
                fprint(structure_file, '\t');
                fprint_id_noabundance(structure_file, data.info(amp_id), parameters.opt_usearch_abundance);
                fprint(structure_file, '\t');
                fprint_integer(structure_file, link_differences);
                fprint(structure_file, '\t');
                fprint_integer(structure_file, cluster_no + 1);
                fprint(structure_file, '\t');
                fprint_integer(structure_file, ampinfo_v[amp_id].generation);
                fprint(structure_file, '\n');
              }
          }

        ++cluster_no;
        progress.update(swarmid + 1);
      }
    progress.done();
  }


  auto write_stats_file(struct Parameters const & parameters,
                        Data const & data,
                        std::vector<struct swarminfo_s> const & swarminfo_v) -> void {
    Progress progress("Writing stats:    ", swarminfo_v.size(), parameters);
    auto * const stats_file = parameters.statsfile.get();

    for (auto const & swarm_info : swarminfo_v) {
      assert(not swarm_info.attached);
      if (swarm_info.attached) {
        continue;
      }
      fprint_integer(stats_file, swarm_info.size);
      fprint(stats_file, '\t');
      fprint_integer(stats_file, swarm_info.mass);
      fprint(stats_file, '\t');
      fprint_id_noabundance(stats_file, data.info(swarm_info.seed), parameters.opt_usearch_abundance);
      fprint(stats_file, '\t');
      fprint_integer(stats_file, data.abundance(swarm_info.seed));
      fprint(stats_file, '\t');
      fprint_integer(stats_file, swarm_info.singletons);
      fprint(stats_file, '\t');
      fprint_integer(stats_file, swarm_info.maxgen);
      fprint(stats_file, '\t');
      fprint_integer(stats_file, swarm_info.maxgen);
      fprint(stats_file, '\n');
      progress.increment();
    }
    progress.done();
  }

} // namespace


auto output_results(struct Parameters const & parameters,
                    Data const & data,
                    std::vector<struct ampinfo_s> const & ampinfo_v,
                    std::vector<struct swarminfo_s> const & swarminfo_v,
                    Overall_stats const & overall_stats) -> void {
  /* dump swarms */
  if (parameters.opt_mothur) {
    write_swarms_mothur_format(parameters, data, ampinfo_v, swarminfo_v, overall_stats);
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
