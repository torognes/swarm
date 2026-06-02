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

#include "algod1_output.h"
#include "../swarm.h"
#include "../db.h"
#include "algod1_internal.h"
#include "nw_aligner.h"
#include "progress.h"
#include "span.h"
#include <algorithm>  // std::sort()
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdio>  // fputc(), fprintf()
#include <iterator>  // std::next()
#include <numeric>  // std::iota
#include <vector>


auto write_network_file(const unsigned int number_of_networks,
                        struct Parameters const & parameters,
                        Data const & data,
                        std::vector<struct ampinfo_s> const & ampinfo_v,
                        std::vector<unsigned int> & network_v) -> void {
  // a network is a cluster with at least two sequences (no singletons)
  Progress progress("Dumping network:  ", number_of_networks, parameters);

  assert(ampinfo_v.size() == data.sequence_count());
  auto counter = 0ULL;
  for (auto const& amplicon: ampinfo_v) {
    const auto link_start = amplicon.link_start;
    const auto link_count = amplicon.link_count;

    auto const neighbours = Span<unsigned int>{std::next(network_v.data(), link_start), link_count};

    // amplicon indexes are already sorted by decreasing abundance
    // then by header in db.cc, so a natural ascending sort here
    // emits neighbours in that ranking order. Earlier dereplication
    // guarantees indexes are distinct.
    std::sort(neighbours.begin(), neighbours.end());

    for (auto const neighbour : neighbours)
      {
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


namespace {

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
                                  std::vector<struct swarminfo_s> const & swarminfo_v,
                                  Overall_stats const & overall_stats) -> void {
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

      auto const seed_seq = data.sequence_view(seed);
      for (auto amp_id = seed_info.next; amp_id != no_swarm; amp_id = ampinfo_v[amp_id].next)
        {
          auto const amp_seq = data.sequence_view(amp_id);

          auto const result = aligner.align(amp_seq, seed_seq);

          std::fprintf(parameters.uclustfile.get(),
                       "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                       cluster_no,
                       amp_seq.length,
                       result.percent_id,
                       result.differences > 0 ? result.cigar_string : "=");

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
        progress.update(swarmid + 1);
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
