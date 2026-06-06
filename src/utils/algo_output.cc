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

#include "algo_output.h"
#include "algo_internal.h"
#include "../swarm.h"
#include "../db.h"
#include "nw_aligner.h"
#include "progress.h"
#include <algorithm>  // std::sort
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc(), fflush
#include <iterator> // std::next
#include <limits>
#include <vector>


namespace {

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

} // namespace


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


auto write_representative_sequences(uint64_t const amplicons,
                                    struct Parameters const & parameters,
                                    Data const & data,
                                    std::vector<struct ampliconinfo_s> const & amps_v) -> void {
  auto seeds = collect_seeds(parameters, data, amplicons, amps_v);
  sort_seeds(parameters, data, seeds);
  write_seeds(parameters, data, seeds);
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
