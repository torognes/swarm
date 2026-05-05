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
#include "qgram.h"
#include "nw.h"
#include "scan.h"
#include "utils/cigar.h"
#include "utils/make_unique.h"
#include "utils/qgram_threadinfo.h"
#include "utils/progress.h"
#include "utils/search_data.h"
#include "utils/score_matrix.h"
#include <algorithm>  // std::min(), std::reverse(), std::for_each
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

  struct ampliconinfo_s {
    unsigned int ampliconid {0};
    unsigned int diffestimate {0}; /* lower bound estimate of dist from initial seed */
    unsigned int swarmid {0};
    unsigned int generation {0};
    unsigned int radius {0}; /* actual diff from initial seed */
  };

  uint64_t swarmed = 0;  // refactoring: reduce scope to algo()?

  struct swarminfo_t {
    uint64_t mass {0};
    unsigned int seed {0};
    int dummy {0}; /* alignment padding only */
  };


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
                     const uint64_t amplicons,
                     std::vector<struct ampliconinfo_s> const & amps_v) -> std::vector<struct swarminfo_t> {
    Progress progress("Collecting seeds:    ", amplicons, parameters);
    assert(swarmed == amplicons);
    std::vector<struct swarminfo_t> seeds(swarmed);  // swarmed == amplicons! Discard swarmed?
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

    // refactoring: replace with std::rotate (rorate rigth by one unit) (reverted: 2x slower)
    // compute reverse iterators starting from the vector's end
    // auto const ssize = static_cast<ptrdiff_t>(amplicons.size());
    // auto const sposition = static_cast<ptrdiff_t>(position);
    // auto const starget = static_cast<ptrdiff_t>(target);
    // auto const first = std::next(amplicons.rbegin(),
    //                              ssize - starget - 1);  // - 1 to include target
    // auto const middle = std::next(first);
    // auto const last = std::next(amplicons.rbegin(), ssize - sposition);
    // std::rotate(first, middle, last);
  }


  auto write_representative_sequences(const uint64_t amplicons,
                                      struct Parameters const & parameters,
                                      Data const & data,
                                      std::vector<struct ampliconinfo_s> const & amps_v) -> void {
    auto seeds = collect_seeds(parameters, data, amplicons, amps_v);
    sort_seeds(parameters, data, seeds);
    write_seeds(parameters, data, seeds);
  }


  auto write_swarms_default_format(const uint64_t amplicons,
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
        const int64_t current_id = amps_v[i].swarmid;
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


  auto write_swarms_mothur_format(const uint64_t amplicons,
                                  const unsigned int swarmid,
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
        const int64_t current_id = amps_v[i].swarmid;
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
} // namespace


auto algo_run(struct Parameters const & parameters,
              Data const & data) -> void {
  auto const score_matrix_63 = create_score_matrix<int64_t>(parameters.penalty_mismatch);

  std::vector<struct Search_data> search_data_v(static_cast<uint64_t>(parameters.opt_threads));
  struct Search_state search_state;
  search_begin(parameters, data, search_state, search_data_v);
  /* start threads */
  assert(parameters.opt_threads <= std::numeric_limits<int>::max());
  auto const search_threads = utils::make_unique<ThreadRunner>(
      static_cast<int>(parameters.opt_threads),
      [&parameters, &data, &search_state](int64_t thread_id) -> void {
        search_worker_core(parameters, data, thread_id, search_state);
      });

  uint64_t largestswarm {0};
  uint64_t maxgenerations {0};

  auto const amplicons = data.sequence_count();
  const uint64_t longestamplicon = data.longest_sequence();

  auto const qgram_store = build_qgram_store(parameters, data);

  std::vector<struct thread_info_s> thread_info_v;
  qgram_diff_init(parameters, qgram_store, thread_info_v);

  std::vector<struct ampliconinfo_s> amps_v(amplicons);
  std::vector<uint64_t> targetampliconids(amplicons);
  std::vector<uint64_t> targetindices(amplicons);
  std::vector<uint64_t> scores_v(amplicons);
  std::vector<uint64_t> diffs_v(amplicons);
  std::vector<uint64_t> alignlengths(amplicons);
  std::vector<uint64_t> qgramamps_v;
  std::vector<uint64_t> qgramdiffs_v(amplicons);
  std::vector<uint64_t> qgramindices_v(amplicons);
  std::vector<uint64_t> hits(amplicons);
  std::vector<unsigned char> directions;
  std::vector<uint64_t> hearray;
  std::vector<char> raw_alignment;
  std::string cigar_string;
  qgramamps_v.reserve(amplicons);
  raw_alignment.reserve(2 * longestamplicon);
  cigar_string.reserve(2 * longestamplicon);

  if (parameters.uclustfile.get() != nullptr)
    {
      directions.resize(longestamplicon * longestamplicon);
      hearray.resize(2 * longestamplicon);
    }

  set_amplicon_ids(amps_v);
  auto const bits = set_bit_mode(parameters);

  uint64_t seeded = 0;

  auto swarmid = 0U;

  Progress progress("Clustering:       ", amplicons, parameters);
  while (seeded < amplicons) {

      /* process each initial seed */

      ++swarmid;
      qgramamps_v.clear();

      uint64_t swarmsize {1};  // a cluster cannot be empty
      uint64_t amplicons_copies {0};  // total abundance of the cluster
      uint64_t singletons {0};
      uint64_t hitcount {0};
      uint64_t maxradius {0};
      uint64_t maxgen {1};  // a cluster can't contain less than 1 generation
      uint64_t seedindex {0};

      seedindex = seeded;
      ++seeded;

      amps_v[seedindex].swarmid = swarmid;

      const uint64_t seedampliconid = amps_v[seedindex].ampliconid;
      hits[hitcount] = seedampliconid;
      ++hitcount;

      auto abundance = data.abundance(seedampliconid);
      amplicons_copies += abundance;
      if (abundance == 1) {
        ++singletons;
      }

      ++swarmed;


      /* find diff estimates between seed and each amplicon in pool */

      uint64_t targetcount = 0;

      // set_list_of_remaining_amplicons
      std::for_each(std::next(amps_v.cbegin(), static_cast<long int>(swarmed)), amps_v.cend(),
          [&parameters, &data, abundance, &qgramamps_v](
              struct ampliconinfo_s const & amplicon) -> void {
            auto const ampliconid = amplicon.ampliconid;
            if ((parameters.opt_no_cluster_breaking) or
                (data.abundance(ampliconid) <= abundance)) {
              qgramamps_v.push_back(ampliconid);
            }
          });
      uint64_t const listlen = qgramamps_v.size();  // temporary refactoring

      qgram_diff_fast(parameters, qgram_store, seedampliconid, listlen, qgramamps_v.data(), qgramdiffs_v.data(), thread_info_v);


      for (auto i = 0ULL; i < listlen; ++i) {
          auto const poolampliconid = qgramamps_v[i];
          auto const diff = qgramdiffs_v[i];
          assert(diff <= std::numeric_limits<unsigned int>::max());
          amps_v[swarmed + i].diffestimate = static_cast<unsigned int>(diff);
          if (diff <= static_cast<uint64_t>(parameters.opt_differences)) {
              targetindices[targetcount] = swarmed + i;
              targetampliconids[targetcount] = poolampliconid;
              ++targetcount;
            }
        }

      if (targetcount > 0) {
          search_do(parameters, data, search_state, seedampliconid, targetcount, targetampliconids.data(),
                    scores_v.data(), diffs_v.data(), alignlengths.data(), bits, search_threads.get());

          for (auto target_id = 0ULL; target_id < targetcount; ++target_id) {
              auto const diff = diffs_v[target_id];

              if (diff <= static_cast<uint64_t>(parameters.opt_differences)) {
                  auto const target = targetindices[target_id];

                  /* move the 'target' to the position ('swarmed')
                     of the first unswarmed amplicon in the pool */
                  move_target_to_first_unswarmed_position(swarmed, target, amps_v);

                  amps_v[swarmed].swarmid = swarmid;
                  amps_v[swarmed].generation = 1;
                  assert(diff <= std::numeric_limits<unsigned int>::max());
                  amps_v[swarmed].radius = static_cast<unsigned int>(diff);
                  maxradius = std::max(diff, maxradius);

                  auto const poolampliconid = amps_v[swarmed].ampliconid;
                  hits[hitcount] = poolampliconid;
                  ++hitcount;

                  if (not parameters.opt_internal_structure.empty()) {
                      data.fprint_id_noabundance(parameters.internal_structure_file.get(),
                                            seedampliconid, parameters.opt_usearch_abundance);
                      std::fprintf(parameters.internal_structure_file.get(), "\t");
                      data.fprint_id_noabundance(parameters.internal_structure_file.get(),
                                            poolampliconid, parameters.opt_usearch_abundance);
                      std::fprintf(parameters.internal_structure_file.get(), "\t%" PRIu64, diff);
                      std::fprintf(parameters.internal_structure_file.get(),
                              "\t%u\t1",
                              swarmid);
                      std::fprintf(parameters.internal_structure_file.get(), "\n");
                    }

                  abundance = data.abundance(poolampliconid);
                  amplicons_copies += abundance;
                  if (abundance == 1) {
                    ++singletons;
                  }

                  ++swarmsize;

                  ++swarmed;
                }
            }


          while (seeded < swarmed) {

              /* process each subseed */

              auto const & subseed = amps_v[seeded];

              ++seeded;

              targetcount = 0;

              auto const subseedabundance = data.abundance(subseed.ampliconid);
              uint64_t subseedlistlen {0};
              for (auto i = swarmed; i < amplicons; ++i) {
                  const uint64_t targetampliconid = amps_v[i].ampliconid;
                  if ((amps_v[i].diffestimate <=
                       subseed.radius + parameters.opt_differences) and
                      ((parameters.opt_no_cluster_breaking) or
                       (data.abundance(targetampliconid)
                        <= subseedabundance))) {
                      qgramamps_v[subseedlistlen] = targetampliconid;
                      qgramindices_v[subseedlistlen] = i;
                      ++subseedlistlen;
                    }
                }

              qgram_diff_fast(parameters, qgram_store, subseed.ampliconid, subseedlistlen, qgramamps_v.data(),
                              qgramdiffs_v.data(), thread_info_v);

              for (auto i = 0ULL; i < subseedlistlen; ++i) {
                if (qgramdiffs_v[i] <= static_cast<uint64_t>(parameters.opt_differences)) {
                  targetindices[targetcount] = qgramindices_v[i];
                  targetampliconids[targetcount] = qgramamps_v[i];
                  ++targetcount;
                }
              }

              if (targetcount == 0) { continue; }

              search_do(parameters, data, search_state, subseed.ampliconid, targetcount, targetampliconids.data(),
                        scores_v.data(), diffs_v.data(), alignlengths.data(), bits, search_threads.get());

              for (auto target_id = 0ULL; target_id < targetcount; ++target_id) {
                  auto const diff = diffs_v[target_id];

                  if (diff > static_cast<uint64_t>(parameters.opt_differences)) { continue; }
                  auto const target = targetindices[target_id];

                  /* find correct position in list */

                  auto const pos = find_correct_position_in_list(swarmed, target, seeded,
                                                                 subseed, amps_v);

                  move_target_to_first_unswarmed_position(pos, target, amps_v);

                  amps_v[pos].swarmid = swarmid;
                  assert(subseed.generation <= std::numeric_limits<unsigned int>::max() - 1);
                  amps_v[pos].generation = subseed.generation + 1;
                  maxgen = std::max<uint64_t>(maxgen, amps_v[pos].generation);
                  assert(subseed.radius <= std::numeric_limits<unsigned int>::max() - diff);
                  amps_v[pos].radius =
                    static_cast<unsigned int>(subseed.radius + diff);
                  maxradius = std::max<uint64_t>(amps_v[pos].radius, maxradius);

                  auto const poolampliconid = amps_v[pos].ampliconid;
                  hits[hitcount] = poolampliconid;
                  ++hitcount;

                  if (not parameters.opt_internal_structure.empty()) {
                    data.fprint_id_noabundance(parameters.internal_structure_file.get(),
                                               subseed.ampliconid,
                                               parameters.opt_usearch_abundance);
                    std::fprintf(parameters.internal_structure_file.get(), "\t");
                    data.fprint_id_noabundance(parameters.internal_structure_file.get(),
                                               poolampliconid,
                                               parameters.opt_usearch_abundance);
                    std::fprintf(parameters.internal_structure_file.get(), "\t%" PRIu64, diff);
                    std::fprintf(parameters.internal_structure_file.get(),
                                 "\t%u\t%u\n",
                                 swarmid, subseed.generation + 1);
                  }

                  abundance = data.abundance(poolampliconid);
                  amplicons_copies += abundance;
                  if (abundance == 1) {
                    ++singletons;
                  }

                  ++swarmsize;

                  ++swarmed;
                }
            }
        }

      largestswarm = std::max(swarmsize, largestswarm);
      maxgenerations = std::max(maxgen, maxgenerations);

      if (parameters.uclustfile.get() != nullptr)
        {
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

            uint64_t nwdiff {0};

            nw(hit_seq.encoded.data(), hit_seq.length, seed_seq.encoded.data(), seed_seq.length,
               score_matrix_63, static_cast<unsigned long int>(parameters.penalty_gapopen),
               static_cast<unsigned long int>(parameters.penalty_gapextend),
               nwdiff, directions, hearray, raw_alignment);

            // backtracking produces a reversed alignment (starting from the end)
            std::reverse(raw_alignment.begin(), raw_alignment.end());
            compress_alignment_to_cigar(raw_alignment, cigar_string);

            // loosing precision when converting raw_alignment.size() and
            // nwdiff to double is not an issue, no need to add assertions
            auto const nwalignmentlength = static_cast<double>(raw_alignment.size());
            auto const differences = static_cast<double>(nwdiff);
            auto const percentid = 100.0 * (nwalignmentlength - differences) / nwalignmentlength;

            std::fprintf(parameters.uclustfile.get(), "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                         swarmid - 1, hit_seq.length, percentid,
                         nwdiff > 0 ? cigar_string.data() : "=");

            data.fprint_id(parameters.uclustfile.get(), hit, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
            std::fprintf(parameters.uclustfile.get(), "\t");
            data.fprint_id(parameters.uclustfile.get(), seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
            std::fprintf(parameters.uclustfile.get(), "\n");
            std::fflush(parameters.uclustfile.get());

            raw_alignment.clear();
            cigar_string.clear();
          }
        }


      if (parameters.statsfile.get() != nullptr) {
        abundance = data.abundance(seedampliconid);

        std::fprintf(parameters.statsfile.get(), "%" PRIu64 "\t%" PRIu64 "\t",
                     swarmsize, amplicons_copies);
        data.fprint_id_noabundance(parameters.statsfile.get(), seedampliconid, parameters.opt_usearch_abundance);
        std::fprintf(parameters.statsfile.get(),
                     "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                     abundance, singletons, maxgen, maxradius);
      }
      progress.update(seeded);
  }
  progress.done();

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

  qgram_diff_done();

  std::fprintf(parameters.logfile, "\n");

  std::fprintf(parameters.logfile, "Number of swarms:  %u\n", swarmid);

  std::fprintf(parameters.logfile, "Largest swarm:     %" PRIu64 "\n", largestswarm);

  std::fprintf(parameters.logfile, "Max generations:   %" PRIu64 "\n", maxgenerations);

  search_end(search_state);
}

