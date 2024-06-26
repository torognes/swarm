/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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
#include "utils/qgram_threadinfo.h"
#include "utils/progress.h"
#include "utils/search_data.h"
#include "utils/seqinfo.h"
#include "utils/score_matrix.h"
#include <algorithm>  // std::min(), std::reverse()
#include <cassert>
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc(), fflush
#include <cstdlib>  // qsort()
#include <cstring>  // strcmp
#include <iterator> // next
#include <limits>
#include <memory>  // unique pointer
#include <string>
#include <vector>

#ifndef PRIu64
#ifdef _WIN32
#define PRIu64 "I64u"
#else
constexpr char PRIu64[] = "lu";
#endif
#endif

#ifndef PRId64
#ifdef _WIN32
#define PRId64 "I64d"
#else
constexpr char PRId64[] = "ld";
#endif
#endif


static uint64_t count_comparisons_8;
static uint64_t count_comparisons_16;

static uint64_t targetcount;

struct ampliconinfo_s
{
  unsigned int ampliconid;
  unsigned int diffestimate; /* lower bound estimate of dist from initial seed */
  unsigned int swarmid;
  unsigned int generation;
  unsigned int radius; /* actual diff from initial seed */
};

static uint64_t swarmed;
static uint64_t seeded;

struct swarminfo_t
{
  uint64_t mass {0};
  unsigned int seed {0};
  int dummy {0}; /* alignment padding only */
};


auto collect_seeds(const uint64_t amplicons,
                   std::vector<struct ampliconinfo_s> & amps_v) -> std::vector<struct swarminfo_t> {
  progress_init("Collecting seeds:    ", amplicons);
  std::vector<struct swarminfo_t> seeds(swarmed);  // swarmed == amplicons! Discard swarmed?
  auto swarmcount = 0UL;
  uint64_t mass = 0;
  auto previous_id = amps_v[0].swarmid;
  auto seed = amps_v[0].ampliconid;
  mass += db_getabundance(seed);
  for(auto i = 1ULL; i < amplicons; ++i)
    {
      const auto current_id = amps_v[i].swarmid;
      if (current_id != previous_id)
        {
          seeds[swarmcount].seed = seed;  // update previous
          seeds[swarmcount].mass = mass;
          ++swarmcount;
          mass = 0;
          seed = amps_v[i].ampliconid;
        }
      mass += db_getabundance(amps_v[i].ampliconid);
      previous_id = current_id;
      progress_update(i);
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
                std::vector<struct swarminfo_t> & seeds) -> void {
  progress_init("Sorting seeds:    ", seeds.size());

  auto compare_seeds = [](struct swarminfo_t const& lhs,
                          struct swarminfo_t const& rhs) -> bool {
    // sort by decreasing mass...
    if (lhs.mass > rhs.mass) {
      return true;
    }
    if (lhs.mass < rhs.mass) {
      return false;
    }
    // ...then ties are sorted by label (alphabetical order)
    auto * const lhs_header = db_getheader(lhs.seed);
    auto * const rhs_header = db_getheader(rhs.seed);
    const auto results = std::strcmp(lhs_header, rhs_header);
    return results == -1;
  };

  std::sort(seeds.begin(), seeds.end(), compare_seeds);
  progress_done(parameters);
}


auto write_seeds(struct Parameters const & parameters,
                 std::vector<struct swarminfo_t> const & seeds) -> void {
  progress_init("Writing seeds:    ", seeds.size());
  auto ticker = 0ULL;  // refactoring: C++20 move ticker to range-loop init-statement
  for(const auto& seed: seeds) {
      const auto swarm_mass = seed.mass;
      const auto swarm_seed = seed.seed;

      std::fprintf(parameters.seeds_file, ">");
      fprint_id_with_new_abundance(parameters.seeds_file, swarm_seed, swarm_mass, parameters.opt_usearch_abundance);
      std::fprintf(parameters.seeds_file, "\n");
      db_fprintseq(parameters.seeds_file, swarm_seed);
      progress_update(ticker);
      ++ticker;
  }
  progress_done(parameters);
}


auto write_representative_sequences(const uint64_t amplicons,
                                    struct Parameters const & parameters,
                                    std::vector<struct ampliconinfo_s> & amps_v) -> void {
  auto seeds = collect_seeds(amplicons, amps_v);
  sort_seeds(parameters, seeds);
  write_seeds(parameters, seeds);
}


auto write_swarms_default_format(const uint64_t amplicons,
                                 struct Parameters const & parameters,
                                 std::vector<struct ampliconinfo_s> & amps_v) -> void {
  /* native swarm output */
  static constexpr char sepchar {' '};  /* usually a space */
  static constexpr char sep_swarms {'\n'};

  fprint_id(parameters.outfile, amps_v[0].ampliconid,
            parameters.opt_usearch_abundance, parameters.opt_append_abundance);
  int64_t previous_id = amps_v[0].swarmid;

  for(auto i = 1ULL; i < amplicons; ++i)
    {
      const int64_t current_id = amps_v[i].swarmid;
      if (current_id == previous_id) {
        std::fputc(sepchar, parameters.outfile);
      }
      else {
        std::fputc(sep_swarms, parameters.outfile);
      }
      fprint_id(parameters.outfile, amps_v[i].ampliconid,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      previous_id = current_id;
    }
  std::fputc('\n', parameters.outfile);
}


auto write_swarms_mothur_format(const uint64_t amplicons,
                                const unsigned int swarmid,
                                struct Parameters const & parameters,
                                std::vector<struct ampliconinfo_s> & amps_v) -> void {
  /* mothur list file output */
  static constexpr char sep_amplicons {','};
  static constexpr char sep_swarms {'\t'};

  std::fprintf(parameters.outfile, "swarm_%" PRId64 "\t%u\t", parameters.opt_differences, swarmid);

  fprint_id(parameters.outfile, amps_v[0].ampliconid,
            parameters.opt_usearch_abundance, parameters.opt_append_abundance);
  int64_t previous_id = amps_v[0].swarmid;

  for(auto i = 1ULL; i < amplicons; ++i)
    {
      const int64_t current_id = amps_v[i].swarmid;
      if (current_id == previous_id) {
        std::fputc(sep_amplicons, parameters.outfile);
      }
      else {
        std::fputc(sep_swarms, parameters.outfile);
      }
      fprint_id(parameters.outfile, amps_v[i].ampliconid,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      previous_id = current_id;
    }

  std::fputc('\n', parameters.outfile);
}


auto algo_run(struct Parameters const & parameters,
              std::vector<struct seqinfo_s> & seqindex_v) -> void
{
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
  const auto score_matrix_63 = create_score_matrix<int64_t>(parameters.penalty_mismatch);

  std::vector<struct Search_data> search_data_v(static_cast<uint64_t>(parameters.opt_threads));
  search_begin(search_data_v);
  /* start threads */
  assert(parameters.opt_threads <= std::numeric_limits<int>::max());
  const std::unique_ptr<ThreadRunner> search_threads (new ThreadRunner(static_cast<int>(parameters.opt_threads), search_worker_core));

  count_comparisons_8 = 0;
  count_comparisons_16 = 0;

  uint64_t largestswarm {0};
  uint64_t maxgenerations {0};

  const auto amplicons = db_getsequencecount();
  const uint64_t longestamplicon = db_getlongestsequence();

  db_qgrams_init(parameters, seqindex_v);

  std::vector<struct thread_info_s> thread_info_v;
  qgram_diff_init(thread_info_v);

  std::vector<struct ampliconinfo_s> amps_v(amplicons);
  std::vector<uint64_t> targetampliconids(amplicons);
  std::vector<uint64_t> targetindices(amplicons);
  std::vector<uint64_t> scores_v(amplicons);
  std::vector<uint64_t> diffs_v(amplicons);
  std::vector<uint64_t> alignlengths(amplicons);
  std::vector<uint64_t> qgramamps_v(amplicons);
  std::vector<uint64_t> qgramdiffs_v(amplicons);
  std::vector<uint64_t> qgramindices_v(amplicons);
  std::vector<uint64_t> hits(amplicons);

  auto const diff_saturation
    = static_cast<uint64_t>(std::min(uint8_max / parameters.penalty_mismatch,
                                     uint8_max / (parameters.penalty_gapopen +
                                                  parameters.penalty_gapextend)));

  std::vector<unsigned char> directions;
  std::vector<uint64_t> hearray;
  std::vector<char> raw_alignment;
  std::string cigar_string;
  raw_alignment.reserve(2 * longestamplicon);
  cigar_string.reserve(2 * longestamplicon);

  if (parameters.uclustfile != nullptr)
    {
      directions.resize(longestamplicon * longestamplicon);
      hearray.resize(2 * longestamplicon);
    }

  /* set ampliconid for all */
  for(auto i = 0U; i < amplicons; ++i) {
    amps_v[i].ampliconid = i;
  }

  /* always search in 8 bit mode unless resolution is very high */
  static constexpr auto bit_mode_8 = 8;
  static constexpr auto bit_mode_16 = 16;
  auto bits = bit_mode_8;

  if (static_cast<uint64_t>(parameters.opt_differences) > diff_saturation) {
    bits = bit_mode_16;
  }

#ifdef __aarch64__
#if !defined(DEBUG) || !defined(COVERAGE)
  /* always use 16-bit version on aarch64 because it is faster */
  bits = bit_mode_16;
#endif
#endif

  seeded = 0;
  swarmed = 0;

  auto swarmid = 0U;

  progress_init("Clustering:       ", amplicons);
  while (seeded < amplicons)
    {

      /* process each initial seed */

      ++swarmid;

      uint64_t swarmsize {0};
      uint64_t amplicons_copies {0};
      uint64_t singletons {0};
      uint64_t hitcount {0};
      uint64_t maxradius {0};
      uint64_t maxgen {1};  // a cluster can't contain less than 1 generation
      uint64_t seedindex {0};

      seedindex = seeded;
      ++seeded;

      amps_v[seedindex].swarmid = swarmid;
      amps_v[seedindex].generation = 0;
      amps_v[seedindex].radius = 0;

      const uint64_t seedampliconid = amps_v[seedindex].ampliconid;
      hits[hitcount] = seedampliconid;
      ++hitcount;

      auto abundance = db_getabundance(seedampliconid);
      amplicons_copies += abundance;
      if (abundance == 1) {
        ++singletons;
      }

      swarmsize = 1;
      ++swarmed;


      /* find diff estimates between seed and each amplicon in pool */

      targetcount = 0;

      uint64_t listlen {0};

      for(auto i = 0ULL; i < amplicons - swarmed; ++i)
        {
          const auto ampid = amps_v[swarmed + i].ampliconid;
          if ((parameters.opt_no_cluster_breaking) or (db_getabundance(ampid) <= abundance))
            {
              qgramamps_v[i] = ampid;
              ++listlen;
            }
        }

      qgram_diff_fast(seedampliconid, listlen, qgramamps_v.data(), qgramdiffs_v.data(), thread_info_v);


      for(auto i = 0ULL; i < listlen; ++i)
        {
          const auto poolampliconid = qgramamps_v[i];
          const auto diff = qgramdiffs_v[i];
          assert(diff <= std::numeric_limits<unsigned int>::max());
          amps_v[swarmed + i].diffestimate = static_cast<unsigned int>(diff);
          if (diff <= static_cast<uint64_t>(parameters.opt_differences))
            {
              targetindices[targetcount] = swarmed + i;
              targetampliconids[targetcount] = poolampliconid;
              ++targetcount;
            }
        }

      if (targetcount > 0)
        {
          search_do(seedampliconid, targetcount, targetampliconids.data(),
                    scores_v.data(), diffs_v.data(), alignlengths.data(), bits, search_threads.get());

          if (bits == bit_mode_8) {
            count_comparisons_8 += targetcount;
          }
          else {
            count_comparisons_16 += targetcount;
          }

          for(auto target_id = 0ULL; target_id < targetcount; ++target_id)
            {
              const auto diff = diffs_v[target_id];

              if (diff <= static_cast<uint64_t>(parameters.opt_differences))
                {
                  const auto target = targetindices[target_id];

                  /* move the 'target' to the position ('swarmed')
                     of the first unswarmed amplicon in the pool */

                  if (swarmed < target)
                    {
                      auto const temp = amps_v[target];  // refactoring: static?
                      for(auto j = target; j > swarmed; --j)
                        {
                          amps_v[j] = amps_v[j - 1];
                        }
                      amps_v[swarmed] = temp;
                    }

                  amps_v[swarmed].swarmid = swarmid;
                  amps_v[swarmed].generation = 1;
                  assert(diff <= std::numeric_limits<unsigned int>::max());
                  amps_v[swarmed].radius = static_cast<unsigned int>(diff);
                  maxradius = std::max(diff, maxradius);

                  const auto poolampliconid = amps_v[swarmed].ampliconid;
                  hits[hitcount] = poolampliconid;
                  ++hitcount;

                  if (not parameters.opt_internal_structure.empty())
                    {
                      fprint_id_noabundance(parameters.internal_structure_file,
                                            seedampliconid, parameters.opt_usearch_abundance);
                      std::fprintf(parameters.internal_structure_file, "\t");
                      fprint_id_noabundance(parameters.internal_structure_file,
                                            poolampliconid, parameters.opt_usearch_abundance);
                      std::fprintf(parameters.internal_structure_file, "\t%" PRIu64, diff);
                      std::fprintf(parameters.internal_structure_file,
                              "\t%u\t1",
                              swarmid);
                      std::fprintf(parameters.internal_structure_file, "\n");
                    }

                  abundance = db_getabundance(poolampliconid);
                  amplicons_copies += abundance;
                  if (abundance == 1) {
                    ++singletons;
                  }

                  ++swarmsize;

                  ++swarmed;
                }
            }


          while (seeded < swarmed)
            {

              /* process each subseed */

              auto subseedampliconid = 0U;
              auto subseedradius = 0U;

              uint64_t subseedindex {0};
              uint64_t subseedgeneration {0};
              uint64_t subseedabundance {0};

              subseedindex = seeded;
              subseedampliconid = amps_v[subseedindex].ampliconid;
              subseedradius = amps_v[subseedindex].radius;
              subseedgeneration = amps_v[subseedindex].generation;
              subseedabundance = db_getabundance(subseedampliconid);

              ++seeded;

              targetcount = 0;

              uint64_t subseedlistlen {0};
              for(auto i = swarmed; i < amplicons; ++i)
                {
                  const uint64_t targetampliconid = amps_v[i].ampliconid;
                  if ((amps_v[i].diffestimate <=
                       subseedradius + parameters.opt_differences) and
                      ((parameters.opt_no_cluster_breaking) or
                       (db_getabundance(targetampliconid)
                        <= subseedabundance)))
                    {
                      qgramamps_v[subseedlistlen] = targetampliconid;
                      qgramindices_v[subseedlistlen] = i;
                      ++subseedlistlen;
                    }
                }

              qgram_diff_fast(subseedampliconid, subseedlistlen, qgramamps_v.data(),
                              qgramdiffs_v.data(), thread_info_v);

              for(auto i = 0ULL; i < subseedlistlen; ++i) {
                if (qgramdiffs_v[i] <= static_cast<uint64_t>(parameters.opt_differences))
                  {
                    targetindices[targetcount] = qgramindices_v[i];
                    targetampliconids[targetcount] = qgramamps_v[i];
                    ++targetcount;
                  }
              }

              if (targetcount > 0)
                {
                  search_do(subseedampliconid, targetcount, targetampliconids.data(),
                            scores_v.data(), diffs_v.data(), alignlengths.data(), bits, search_threads.get());

                  if (bits == bit_mode_8) {
                    count_comparisons_8 += targetcount;
                  }
                  else {
                    count_comparisons_16 += targetcount;
                  }

                  for(auto target_id = 0ULL; target_id < targetcount; ++target_id)
                    {
                      const auto diff = diffs_v[target_id];

                      if (diff <= static_cast<uint64_t>(parameters.opt_differences))
                        {
                          const auto target = targetindices[target_id];

                          /* find correct position in list */

                          /* move the 'target' to the position (swarmed)
                             of the first unswarmed amplicon in the pool
                             then move the target further into the swarmed
                             but unseeded part of the list, so that the
                             swarmed amplicons are ordered by id */

                          const uint64_t targetampliconid = amps_v[target].ampliconid;
                          auto pos = swarmed;

                          while ((pos > seeded) and
                                 (amps_v[pos - 1].ampliconid > targetampliconid) and
                                 (amps_v[pos - 1].generation > subseedgeneration)) {
                            --pos;
                          }

                          if (pos < target)
                            {
                              auto const temp = amps_v[target];  // refactoring: static?
                              for(auto j = target; j > pos; --j)
                                {
                                  amps_v[j] = amps_v[j - 1];
                                }
                              amps_v[pos] = temp;
                            }

                          amps_v[pos].swarmid = swarmid;
                          assert(subseedgeneration + 1 <= std::numeric_limits<unsigned int>::max());
                          amps_v[pos].generation =
                              static_cast<unsigned int>(subseedgeneration + 1);
                          maxgen = std::max<uint64_t>(maxgen, amps_v[pos].generation);
                          assert(subseedradius + diff <= std::numeric_limits<unsigned int>::max());
                          amps_v[pos].radius =
                              static_cast<unsigned int>(subseedradius + diff);
                          maxradius = std::max<uint64_t>(amps_v[pos].radius, maxradius);

                          const auto poolampliconid = amps_v[pos].ampliconid;
                          hits[hitcount] = poolampliconid;
                          ++hitcount;

                          if (not parameters.opt_internal_structure.empty())
                            {
                              fprint_id_noabundance(parameters.internal_structure_file,
                                                    subseedampliconid,
                                                    parameters.opt_usearch_abundance);
                              std::fprintf(parameters.internal_structure_file, "\t");
                              fprint_id_noabundance(parameters.internal_structure_file,
                                                    poolampliconid,
                                                    parameters.opt_usearch_abundance);
                              std::fprintf(parameters.internal_structure_file, "\t%" PRIu64, diff);
                              std::fprintf(parameters.internal_structure_file,
                                      "\t%u\t%" PRIu64,
                                      swarmid, subseedgeneration + 1);
                              std::fprintf(parameters.internal_structure_file, "\n");
                            }

                          abundance = db_getabundance(poolampliconid);
                          amplicons_copies += abundance;
                          if (abundance == 1) {
                            ++singletons;
                          }

                          ++swarmsize;

                          ++swarmed;
                        }
                    }
                }
            }
        }

      largestswarm = std::max(swarmsize, largestswarm);
      maxgenerations = std::max(maxgen, maxgenerations);

      if (parameters.uclustfile != nullptr)
        {
          std::fprintf(parameters.uclustfile, "C\t%u\t%" PRIu64 "\t*\t*\t*\t*\t*\t",
                  swarmid-1, swarmsize);
          fprint_id(parameters.uclustfile, seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile, "\t*\n");

          std::fprintf(parameters.uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                  swarmid-1, db_getsequencelen(seedampliconid));
          fprint_id(parameters.uclustfile, seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile, "\t*\n");
          std::fflush(parameters.uclustfile);

          for(auto i = 1ULL; i < hitcount; ++i)
            {
              const auto hit = hits[i];

              auto * dseq = db_getsequence(hit);
              const auto dlen = db_getsequencelen(hit);
              auto * qseq = db_getsequence(seedampliconid);
              const auto qlen = db_getsequencelen(seedampliconid);

              uint64_t nwdiff {0};

              nw(dseq, dlen, qseq, qlen,
                 score_matrix_63, static_cast<unsigned long int>(parameters.penalty_gapopen),
                 static_cast<unsigned long int>(parameters.penalty_gapextend),
                 nwdiff, directions, hearray, raw_alignment);

              // backtracking produces a reversed alignment (starting from the end)
              std::reverse(raw_alignment.begin(), raw_alignment.end());
              compress_alignment_to_cigar(raw_alignment, cigar_string);

              // loosing precision when converting raw_alignment.size() and
              // nwdiff to double is not an issue, no need to add assertions
              const auto nwalignmentlength = static_cast<double>(raw_alignment.size());
              const auto differences = static_cast<double>(nwdiff);
              const double percentid = 100.0 * (nwalignmentlength - differences) / nwalignmentlength;

              std::fprintf(parameters.uclustfile, "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                      swarmid-1, db_getsequencelen(hit), percentid,
                      nwdiff > 0 ? cigar_string.data() : "=");

              fprint_id(parameters.uclustfile, hit, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
              std::fprintf(parameters.uclustfile, "\t");
              fprint_id(parameters.uclustfile, seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
              std::fprintf(parameters.uclustfile, "\n");
              std::fflush(parameters.uclustfile);

              raw_alignment.clear();
              cigar_string.clear();
            }

        }


      if (parameters.statsfile != nullptr)
        {
          abundance = db_getabundance(seedampliconid);

          std::fprintf(parameters.statsfile, "%" PRIu64 "\t%" PRIu64 "\t",
                  swarmsize, amplicons_copies);
          fprint_id_noabundance(parameters.statsfile, seedampliconid, parameters.opt_usearch_abundance);
          std::fprintf(parameters.statsfile,
                  "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                  abundance, singletons, maxgen, maxradius);
        }
      progress_update(seeded);
    }
  progress_done(parameters);

  /* output swarms */
  if (amplicons > 0) {
    if (parameters.opt_mothur) {
      write_swarms_mothur_format(amplicons, swarmid, parameters, amps_v);
    }
    else {
      write_swarms_default_format(amplicons, parameters, amps_v);
    }
  }


  /* dump seeds in fasta format with sum of abundances */
  if ((not parameters.opt_seeds.empty()) and (amplicons > 0)) {
    write_representative_sequences(amplicons, parameters, amps_v);
  }

  db_qgrams_done();

  qgram_diff_done();

  std::fprintf(parameters.logfile, "\n");

  std::fprintf(parameters.logfile, "Number of swarms:  %u\n", swarmid);

  std::fprintf(parameters.logfile, "Largest swarm:     %" PRIu64 "\n", largestswarm);

  std::fprintf(parameters.logfile, "Max generations:   %" PRIu64 "\n", maxgenerations);

  search_end();
}

