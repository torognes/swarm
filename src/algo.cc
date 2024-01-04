/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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
#include "util.h"
#include "utils/score_matrix.h"
#include <algorithm>  // std::min()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fputc()
#include <cstdlib>  // qsort()
#include <cstring>  // strcmp
#include <iterator> // next
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
static uint64_t * scores;
static uint64_t * diffs;
static uint64_t * qgramamps;
static uint64_t * qgramdiffs;
static uint64_t * qgramindices;

static struct ampliconinfo_s
{
  unsigned ampliconid;
  unsigned diffestimate; /* lower bound estimate of dist from initial seed */
  unsigned swarmid;
  unsigned generation;
  unsigned radius; /* actual diff from initial seed */
} * amps;

static uint64_t swarmed;
static uint64_t seeded;

struct swarminfo_t
{
  uint64_t mass {0};
  unsigned int seed {0};
  int dummy {0}; /* alignment padding only */
};


auto collect_seeds(const uint64_t amplicons) -> std::vector<struct swarminfo_t> {
  progress_init("Collecting seeds:    ", amplicons);
  std::vector<struct swarminfo_t> seeds(swarmed);  // swarmed == amplicons! Discard swarmed?
  auto swarmcount = 0UL;
  uint64_t mass = 0;
  unsigned int previd = amps[0].swarmid;
  unsigned int seed = amps[0].ampliconid;
  mass += db_getabundance(seed);
  for(auto i = 1ULL; i < amplicons; i++)
    {
      const unsigned int id = amps[i].swarmid;
      if (id != previd)
        {
          seeds[swarmcount].seed = seed;  // update previous
          seeds[swarmcount].mass = mass;
          ++swarmcount;
          mass = 0;
          seed = amps[i].ampliconid;
        }
      mass += db_getabundance(amps[i].ampliconid);
      previd = id;
      progress_update(i);
    }
  seeds[swarmcount].seed = seed;
  seeds[swarmcount].mass = mass;
  ++swarmcount;

  // free some memory
  seeds.erase(std::next(seeds.begin(), static_cast<long int>(swarmcount)), seeds.end());
  seeds.shrink_to_fit();

  return seeds;
}


auto sort_seeds(std::vector<struct swarminfo_t>& seeds) -> void {
  progress_init("Sorting seeds:    ", seeds.size());

  auto compare_seeds = [](struct swarminfo_t const& lhs, struct swarminfo_t const& rhs) -> bool {
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
  progress_done();
}


auto write_seeds(std::vector<struct swarminfo_t> const & seeds,
                 bool const opt_usearch_abundance) -> void {
  progress_init("Writing seeds:    ", seeds.size());
  auto ticker = 0ULL;  // refactoring: C++20 move ticker to range-loop init-statement
  for(const auto& seed: seeds) {
      const auto swarm_mass = seed.mass;
      const auto swarm_seed = seed.seed;

      std::fprintf(fp_seeds, ">");
      fprint_id_with_new_abundance(fp_seeds, swarm_seed, swarm_mass, opt_usearch_abundance);
      std::fprintf(fp_seeds, "\n");
      db_fprintseq(fp_seeds, swarm_seed);
      progress_update(ticker);
      ++ticker;
  }
  progress_done();
}


auto write_representative_sequences(const uint64_t amplicons,
                                    struct Parameters const & parameters) -> void {
  auto seeds = collect_seeds(amplicons);
  sort_seeds(seeds);
  write_seeds(seeds, parameters.opt_usearch_abundance);
}


auto write_swarms_default_format(const uint64_t amplicons,
                                 struct Parameters const & parameters) -> void {
  /* native swarm output */
  static constexpr char sep_amplicons {sepchar};  /* usually a space */
  static constexpr char sep_swarms {'\n'};

  fprint_id(outfile, amps[0].ampliconid,
            parameters.opt_usearch_abundance, parameters.opt_append_abundance);
  int64_t previd = amps[0].swarmid;

  for(auto i = 1ULL; i < amplicons; i++)
    {
      const int64_t id = amps[i].swarmid;
      if (id == previd) {
        std::fputc(sep_amplicons, outfile);
      }
      else {
        std::fputc(sep_swarms, outfile);
      }
      fprint_id(outfile, amps[i].ampliconid,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      previd = id;
    }
  std::fputc('\n', outfile);
}


auto write_swarms_mothur_format(const uint64_t amplicons,
                                const unsigned int swarmid,
                                struct Parameters const & parameters) -> void {
  /* mothur list file output */
  static constexpr char sep_amplicons {','};
  static constexpr char sep_swarms {'\t'};

  std::fprintf(outfile, "swarm_%" PRId64 "\t%u\t", parameters.opt_differences, swarmid);

  fprint_id(outfile, amps[0].ampliconid,
            parameters.opt_usearch_abundance, parameters.opt_append_abundance);
  int64_t previd = amps[0].swarmid;

  for(auto i = 1ULL; i < amplicons; i++)
    {
      const int64_t id = amps[i].swarmid;
      if (id == previd) {
        std::fputc(sep_amplicons, outfile);
      }
      else {
        std::fputc(sep_swarms, outfile);
      }
      fprint_id(outfile, amps[i].ampliconid,
                parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      previd = id;
    }

  std::fputc('\n', outfile);
}


auto algo_run(struct Parameters const & parameters) -> void
{
  const auto score_matrix_63 = create_score_matrix<int64_t>(parameters.penalty_mismatch);

  search_begin();

  count_comparisons_8 = 0;
  count_comparisons_16 = 0;


  uint64_t largestswarm {0};

  uint64_t maxgenerations {0};

  const uint64_t amplicons = db_getsequencecount();
  const uint64_t longestamplicon = db_getlongestsequence();

  db_qgrams_init();

  qgram_diff_init();

  amps = new struct ampliconinfo_s[amplicons];
  std::vector<uint64_t> targetampliconids(amplicons);
  std::vector<uint64_t> targetindices(amplicons);
  scores = new uint64_t[amplicons];
  diffs = new uint64_t[amplicons];
  std::vector<uint64_t> alignlengths(amplicons);
  qgramamps = new uint64_t[amplicons];
  qgramdiffs = new uint64_t[amplicons];
  qgramindices = new uint64_t[amplicons];
  std::vector<uint64_t> hits(amplicons);

  auto diff_saturation
    = static_cast<uint64_t>(std::min(UINT8_MAX / parameters.penalty_mismatch,
                                     UINT8_MAX / (penalty_gapopen +
                                                  penalty_gapextend)));

  std::vector<unsigned char> dir;
  std::vector<uint64_t> hearray;

  if (uclustfile != nullptr)
    {
      dir.resize(longestamplicon * longestamplicon);
      hearray.resize(2 * longestamplicon);
    }

  /* set ampliconid for all */
  for(auto i = 0U; i < amplicons; i++) {
    amps[i].ampliconid = i;
  }

  /* always search in 8 bit mode unless resolution is very high */
  static constexpr int bit_mode_8 {8};
  static constexpr int bit_mode_16 {16};
  int bits {bit_mode_8};

  if (static_cast<uint64_t>(parameters.opt_differences) > diff_saturation) {
    bits = bit_mode_16;
  }

#ifdef __aarch64__
  /* always use 16-bit version on aarch64 because it is faster */
  bits = bit_mode_16;
#endif

  seeded = 0;
  swarmed = 0;

  unsigned int swarmid {0};

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
      uint64_t maxgen {0};
      uint64_t seedindex {0};

      seedindex = seeded;
      ++seeded;

      amps[seedindex].swarmid = swarmid;
      amps[seedindex].generation = 0;
      amps[seedindex].radius = 0;

      const uint64_t seedampliconid = amps[seedindex].ampliconid;
      hits[hitcount++] = seedampliconid;

      uint64_t abundance = db_getabundance(seedampliconid);
      amplicons_copies += abundance;
      if (abundance == 1) {
        ++singletons;
      }

      swarmsize = 1;
      ++swarmed;


      /* find diff estimates between seed and each amplicon in pool */

      targetcount = 0;

      uint64_t listlen {0};

      for(auto i = 0ULL; i < amplicons - swarmed; i++)
        {
          const unsigned int ampid = amps[swarmed + i].ampliconid;
          if ((opt_no_otu_breaking) or (db_getabundance(ampid) <= abundance))
            {
              qgramamps[i] = ampid;
              ++listlen;
            }
        }

      qgram_diff_fast(seedampliconid, listlen, qgramamps, qgramdiffs);


      for(auto i = 0ULL; i < listlen; i++)
        {
          const uint64_t poolampliconid = qgramamps[i];
          const uint64_t diff = qgramdiffs[i];
          amps[swarmed + i].diffestimate = static_cast<unsigned int>(diff);
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
                    scores, diffs, alignlengths.data(), bits);

          if (bits == bit_mode_8) {
            count_comparisons_8 += targetcount;
          }
          else {
            count_comparisons_16 += targetcount;
          }

          for(auto t = 0ULL; t < targetcount; t++)
            {
              const uint64_t diff = diffs[t];

              if (diff <= static_cast<uint64_t>(parameters.opt_differences))
                {
                  const uint64_t i = targetindices[t];

                  /* move the target (i) to the position (swarmed)
                     of the first unswarmed amplicon in the pool */

                  if (swarmed < i)
                    {
                      const struct ampliconinfo_s temp = amps[i];
                      for(auto j = i; j > swarmed; j--)
                        {
                          amps[j] = amps[j-1];
                        }
                      amps[swarmed] = temp;
                    }

                  amps[swarmed].swarmid = swarmid;
                  amps[swarmed].generation = 1;
                  if (maxgen < 1) {
                    maxgen = 1;
                  }
                  amps[swarmed].radius = static_cast<unsigned int>(diff);
                  if (diff > maxradius) {
                    maxradius = diff;
                  }

                  const unsigned int poolampliconid = amps[swarmed].ampliconid;
                  hits[hitcount++] = poolampliconid;

                  if (not parameters.opt_internal_structure.empty())
                    {
                      fprint_id_noabundance(internal_structure_file,
                                            seedampliconid, parameters.opt_usearch_abundance);
                      std::fprintf(internal_structure_file, "\t");
                      fprint_id_noabundance(internal_structure_file,
                                            poolampliconid, parameters.opt_usearch_abundance);
                      std::fprintf(internal_structure_file, "\t%" PRIu64, diff);
                      std::fprintf(internal_structure_file,
                              "\t%u\t1",
                              swarmid);
                      std::fprintf(internal_structure_file, "\n");
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

              unsigned int subseedampliconid {0};
              unsigned int subseedradius {0};

              uint64_t subseedindex {0};
              uint64_t subseedgeneration {0};
              uint64_t subseedabundance {0};

              subseedindex = seeded;
              subseedampliconid = amps[subseedindex].ampliconid;
              subseedradius = amps[subseedindex].radius;
              subseedgeneration = amps[subseedindex].generation;
              subseedabundance = db_getabundance(subseedampliconid);

              ++seeded;

              targetcount = 0;

              uint64_t subseedlistlen {0};
              for(auto i = swarmed; i < amplicons; i++)
                {
                  const uint64_t targetampliconid = amps[i].ampliconid;
                  if ((amps[i].diffestimate <=
                       subseedradius + parameters.opt_differences) and
                      ((opt_no_otu_breaking) or
                       (db_getabundance(targetampliconid)
                        <= subseedabundance)))
                    {
                      qgramamps[subseedlistlen] = targetampliconid;
                      qgramindices[subseedlistlen] = i;
                      ++subseedlistlen;
                    }
                }

              qgram_diff_fast(subseedampliconid, subseedlistlen, qgramamps,
                              qgramdiffs);

              for(auto i = 0ULL; i < subseedlistlen; i++) {
                if (qgramdiffs[i] <= static_cast<uint64_t>(parameters.opt_differences))
                  {
                    targetindices[targetcount] = qgramindices[i];
                    targetampliconids[targetcount] = qgramamps[i];
                    ++targetcount;
                  }
              }

              if (targetcount > 0)
                {
                  search_do(subseedampliconid, targetcount, targetampliconids.data(),
                            scores, diffs, alignlengths.data(), bits);

                  if (bits == bit_mode_8) {
                    count_comparisons_8 += targetcount;
                  }
                  else {
                    count_comparisons_16 += targetcount;
                  }

                  for(auto t = 0ULL; t < targetcount; t++)
                    {
                      const uint64_t diff = diffs[t];

                      if (diff <= static_cast<uint64_t>(parameters.opt_differences))
                        {
                          const uint64_t i = targetindices[t];

                          /* find correct position in list */

                          /* move the target (i) to the position (swarmed)
                             of the first unswarmed amplicon in the pool
                             then move the target further into the swarmed
                             but unseeded part of the list, so that the
                             swarmed amplicons are ordered by id */

                          const uint64_t targetampliconid = amps[i].ampliconid;
                          uint64_t pos = swarmed;

                          while ((pos > seeded) and
                                 (amps[pos-1].ampliconid > targetampliconid) and
                                 (amps[pos-1].generation > subseedgeneration)) {
                            --pos;
                          }

                          if (pos < i)
                            {
                              const struct ampliconinfo_s temp = amps[i];
                              for(auto j = i; j > pos; j--)
                                {
                                  amps[j] = amps[j-1];
                                }
                              amps[pos] = temp;
                            }

                          amps[pos].swarmid = swarmid;
                          amps[pos].generation
                            = static_cast<unsigned int>(subseedgeneration + 1);
                          if (maxgen < amps[pos].generation) {
                            maxgen = amps[pos].generation;
                          }
                          amps[pos].radius
                            = static_cast<unsigned int>(subseedradius + diff);
                          if (amps[pos].radius > maxradius) {
                            maxradius = amps[pos].radius;
                          }

                          const unsigned int poolampliconid = amps[pos].ampliconid;
                          hits[hitcount++] = poolampliconid;

                          if (not parameters.opt_internal_structure.empty())
                            {
                              fprint_id_noabundance(internal_structure_file,
                                                    subseedampliconid,
                                                    parameters.opt_usearch_abundance);
                              std::fprintf(internal_structure_file, "\t");
                              fprint_id_noabundance(internal_structure_file,
                                                    poolampliconid,
                                                    parameters.opt_usearch_abundance);
                              std::fprintf(internal_structure_file, "\t%" PRIu64, diff);
                              std::fprintf(internal_structure_file,
                                      "\t%u\t%" PRIu64,
                                      swarmid, subseedgeneration + 1);
                              std::fprintf(internal_structure_file, "\n");
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

      if (swarmsize > largestswarm) {
        largestswarm = swarmsize;
      }

      if (maxgen > maxgenerations) {
        maxgenerations = maxgen;
      }


      if (uclustfile != nullptr)
        {
          std::fprintf(uclustfile, "C\t%u\t%" PRIu64 "\t*\t*\t*\t*\t*\t",
                  swarmid-1, swarmsize);
          fprint_id(uclustfile, seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(uclustfile, "\t*\n");

          std::fprintf(uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                  swarmid-1, db_getsequencelen(seedampliconid));
          fprint_id(uclustfile, seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(uclustfile, "\t*\n");
          fflush(uclustfile);

          for(auto i = 1ULL; i < hitcount; i++)
            {
              const auto hit = hits[i];

              char * dseq = db_getsequence(hit);
              const auto dlen = db_getsequencelen(hit);
              char * qseq = db_getsequence(seedampliconid);
              const auto qlen = db_getsequencelen(seedampliconid);

              int64_t nwscore {0};
              int64_t nwdiff {0};
              char * nwalignment {nullptr};
              int64_t nwalignmentlength {0};

              nw(dseq, dlen, qseq, qlen,
                 score_matrix_63, penalty_gapopen, penalty_gapextend,
                 nwscore, nwdiff, nwalignmentlength, & nwalignment,
                 dir, reinterpret_cast<int64_t *>(hearray.data()), 0, 0);

              const double percentid
                = 100.0 * static_cast<double>(nwalignmentlength - nwdiff)
                / static_cast<double>(nwalignmentlength);

              std::fprintf(uclustfile, "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                      swarmid-1, db_getsequencelen(hit), percentid,
                      nwdiff > 0 ? nwalignment : "=");

              fprint_id(uclustfile, hit, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
              std::fprintf(uclustfile, "\t");
              fprint_id(uclustfile, seedampliconid, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
              std::fprintf(uclustfile, "\n");
              fflush(uclustfile);

              if (nwalignment != nullptr) {
                xfree(nwalignment);
              }
              nwalignment = nullptr;
            }

        }


      if (statsfile != nullptr)
        {
          abundance = db_getabundance(seedampliconid);

          std::fprintf(statsfile, "%" PRIu64 "\t%" PRIu64 "\t",
                  swarmsize, amplicons_copies);
          fprint_id_noabundance(statsfile, seedampliconid, parameters.opt_usearch_abundance);
          std::fprintf(statsfile,
                  "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\t%" PRIu64 "\n",
                  abundance, singletons, maxgen, maxradius);
        }
      progress_update(seeded);
    }
  progress_done();

  /* output swarms */
  if (amplicons > 0) {
    if (parameters.opt_mothur) {
      write_swarms_mothur_format(amplicons, swarmid, parameters);
    }
    else {
      write_swarms_default_format(amplicons, parameters);
    }
  }


  /* dump seeds in fasta format with sum of abundances */
  if ((not parameters.opt_seeds.empty()) and (amplicons > 0)) {
    write_representative_sequences(amplicons, parameters);
  }


  delete [] qgramdiffs;
  qgramdiffs = nullptr;
  delete [] qgramamps;
  qgramamps = nullptr;
  delete [] qgramindices;
  qgramindices = nullptr;
  delete [] diffs;
  diffs = nullptr;
  delete [] scores;
  scores = nullptr;
  delete [] amps;
  amps = nullptr;

  db_qgrams_done();

  qgram_diff_done();

  std::fprintf(logfile, "\n");

  std::fprintf(logfile, "Number of swarms:  %u\n", swarmid);

  std::fprintf(logfile, "Largest swarm:     %" PRIu64 "\n", largestswarm);

  std::fprintf(logfile, "Max generations:   %" PRIu64 "\n", maxgenerations);

  search_end();
}

