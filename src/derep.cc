/*
    SWARM

    Copyright (C) 2012-2025 Torbjorn Rognes and Frederic Mahe

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
#include "utils/hashtable_size.h"
#include "utils/nt_codec.h"
#include "utils/progress.h"
#include "zobrist.h"
#include <algorithm>  // sort
#include <cassert>  // assert
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>
#include <cstdio>  // fputc()
#include <cstdlib>  // qsort()
#include <iterator>  // std::next
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


namespace {

  struct bucket
  {
    uint64_t hash = 0;
    unsigned int seqno_first = 0;
    unsigned int seqno_last = 0;
    uint64_t mass = 0;
    unsigned int size = 0;
    unsigned int singletons = 0;
  };

  struct Stats
  {
    int64_t swarmcount = 0;
    uint64_t maxmass = 0;
    unsigned int maxsize = 0U;
  };


  auto sort_seeds(struct Parameters const & parameters,
                  std::vector<struct bucket>& hashtable) -> void {
    progress_init("Sorting:          ", 1);

    auto compare_seeds = [](struct bucket const& lhs,
                            struct bucket const& rhs) -> bool {
      // sort by decreasing mass...
      if (lhs.mass > rhs.mass) {
        return true;
      }
      if (lhs.mass < rhs.mass) {
        return false;
      }
      // ...then ties are sorted by input order
      if (lhs.seqno_first < rhs.seqno_first) {
        return true;
      }
      return false;
    };

    std::sort(hashtable.begin(), hashtable.end(), compare_seeds);
    progress_done(parameters);
  }


  auto release_unused_memory(std::vector<struct bucket>& hashtable,
                             int64_t const swarmcount) -> void {
    hashtable.erase(hashtable.begin() + swarmcount, hashtable.end());
    hashtable.shrink_to_fit();
    assert(hashtable.size() == static_cast<uint64_t>(swarmcount));
  }


  auto write_stats_file(struct Parameters const & parameters,
                        std::vector<struct bucket> const & hashtable) -> void {
    progress_init("Writing stats:    ", hashtable.size());
    auto counter = 0U;
    for(auto const & cluster: hashtable) {
      std::fprintf(parameters.statsfile, "%u\t%" PRIu64 "\t", cluster.size, cluster.mass);
      fprint_id_noabundance(parameters.statsfile, cluster.seqno_first, parameters.opt_usearch_abundance);
      std::fprintf(parameters.statsfile, "\t%" PRIu64 "\t%u\t%u\t%u\n",
                   db_getabundance(cluster.seqno_first),
                   cluster.singletons, 0U, 0U);
      ++counter;
      progress_update(counter);
    }
    progress_done(parameters);
  }


  auto write_structure_file(struct Parameters const & parameters,
                            std::vector<struct bucket> const & hashtable,
                            std::vector<unsigned int> const & nextseqtab) -> void {
    progress_init("Writing structure:", hashtable.size());
    auto counter = 0UL;

    for(auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;
      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          fprint_id_noabundance(parameters.internal_structure_file, seed, parameters.opt_usearch_abundance);
          std::fprintf(parameters.internal_structure_file, "\t");
          fprint_id_noabundance(parameters.internal_structure_file, next_identical, parameters.opt_usearch_abundance);
          std::fprintf(parameters.internal_structure_file, "\t%d\t%lu\t%d\n", 0, counter + 1, 0);
          next_identical = nextseqtab[next_identical];
        }
      ++counter;
      progress_update(counter);
    }
    progress_done(parameters);
  }


  auto write_swarms_uclust_format(struct Parameters const & parameters,
                                  std::vector<struct bucket> const & hashtable,
                                  std::vector<unsigned int> const & nextseqtab) -> void {
    progress_init("Writing UCLUST:   ", hashtable.size());
    auto counter = 0U;

    for(auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;

      std::fprintf(parameters.uclustfile, "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                   counter,
                   cluster.size);
      fprint_id(parameters.uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile, "\t*\n");

      std::fprintf(parameters.uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                   counter,
                   db_getsequencelen(seed));
      fprint_id(parameters.uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile, "\t*\n");

      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          std::fprintf(parameters.uclustfile,
                       "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                       counter,
                       db_getsequencelen(next_identical),
                       100.0,
                       "=");
          fprint_id(parameters.uclustfile, next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile, "\t");
          fprint_id(parameters.uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile, "\n");
          next_identical = nextseqtab[next_identical];
        }
      ++counter;
      progress_update(counter);
    }
    progress_done(parameters);
  }


  auto write_representative_sequences(struct Parameters const & parameters,
                                      std::vector<struct bucket> const & hashtable) -> void {
    progress_init("Writing seeds:    ", hashtable.size());
    auto counter = 0U;
    for(auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;
      std::fprintf(parameters.seeds_file, ">");
      fprint_id_with_new_abundance(parameters.seeds_file, seed, cluster.mass, parameters.opt_usearch_abundance);
      std::fprintf(parameters.seeds_file, "\n");
      db_fprintseq(parameters.seeds_file, seed);
      ++counter;
      progress_update(counter);
    }
    progress_done(parameters);
  }


  auto write_swarms_mothur_format(struct Parameters const & parameters,
                                  std::vector<struct bucket> const & hashtable,
                                  std::vector<unsigned int> const & nextseqtab) -> void {
    progress_init("Writing swarms:   ", hashtable.size());

#ifdef _WIN32
    std::fprintf(parameters.outfile, "swarm_%" PRId64 "\t%llu", parameters.opt_differences, hashtable.size());
#else
    std::fprintf(parameters.outfile, "swarm_%" PRId64 "\t%lu", parameters.opt_differences, hashtable.size());
#endif

    auto counter = 0U;

    for(auto const & cluster: hashtable) {
      // print cluster seed
      auto const seed = cluster.seqno_first;
      std::fputc('\t', parameters.outfile);
      fprint_id(parameters.outfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      // print other cluster members
      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          std::fputc(',', parameters.outfile);
          fprint_id(parameters.outfile, next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          next_identical = nextseqtab[next_identical];
        }

      ++counter;
      progress_update(counter);
    }
    std::fputc('\n', parameters.outfile);

    progress_done(parameters);
  }


  auto write_swarms_default_format(struct Parameters const & parameters,
                                   std::vector<struct bucket> const & hashtable,
                                   std::vector<unsigned int> const & nextseqtab) -> void {
    static constexpr char sepchar {' '};
    progress_init("Writing swarms:   ", hashtable.size());
    auto counter = 0U;

    for(auto const & cluster: hashtable) {
      // print cluster seed
      auto const seed = cluster.seqno_first;
      fprint_id(parameters.outfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      // print other cluster members
      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          std::fputc(sepchar, parameters.outfile);
          fprint_id(parameters.outfile, next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          next_identical = nextseqtab[next_identical];
        }
      std::fputc('\n', parameters.outfile);
      ++counter;
      progress_update(counter);
    }

    progress_done(parameters);
  }


  auto dereplicating(struct Parameters const & parameters,
                     std::vector<struct bucket> & hashtable,
                     std::vector<unsigned int> & nextseqtab)
    -> struct Stats
       {
         progress_init("Dereplicating:    ", nextseqtab.size());

         struct Stats stats;
         const uint64_t derep_hash_mask = hashtable.size() - 1;

         for(auto seqno = 0U; seqno < nextseqtab.size(); ++seqno)
           {
             const unsigned int seqlen = db_getsequencelen(seqno);
             char * seq = db_getsequence(seqno);

             /*
               Find free bucket or bucket for identical sequence.
               Make sure sequences are exactly identical
               in case of any hash collision.
               With 64-bit hashes, there is about 50% chance of a
               collision when the number of sequences is about 5e9.
             */

             const uint64_t hash = zobrist_hash(seq, seqlen);

             uint64_t nth_bucket = hash & derep_hash_mask;
             auto * clusterp = &hashtable[nth_bucket];

             while ((clusterp->mass != 0U) and
                    ((clusterp->hash != hash) or
                     (seqlen != db_getsequencelen(clusterp->seqno_first)) or
                     not std::equal(seq, std::next(seq, nt_bytelength(seqlen)),
                                    db_getsequence(clusterp->seqno_first))
                     )
                    )
               {
                 clusterp = std::next(clusterp);
                 ++nth_bucket;
                 if (nth_bucket >= hashtable.size()) // wrap around the table if we reach the end
                   {
                     nth_bucket = 0;
                     clusterp = hashtable.data();
                   }
               }

             const uint64_t abundance = db_getabundance(seqno);

             if (clusterp->mass != 0U)
               {
                 /* at least one identical sequence already */
                 nextseqtab[clusterp->seqno_last] = seqno;
               }
             else
               {
                 /* no identical sequences yet, start a new cluster */
                 ++stats.swarmcount;
                 clusterp->hash = hash;
                 clusterp->seqno_first = seqno;
                 clusterp->size = 0;
                 clusterp->singletons = 0;
               }

             ++clusterp->size;
             clusterp->seqno_last = seqno;
             clusterp->mass += abundance;

             if (abundance == 1) {
               ++clusterp->singletons;
             }

             stats.maxmass = std::max(clusterp->mass, stats.maxmass);
             stats.maxsize = std::max(clusterp->size, stats.maxsize);

             progress_update(seqno);
           }
         progress_done(parameters);

         return stats;
    }


  auto output_results(struct Parameters const & parameters,
                      std::vector<struct bucket> & hashtable,
                      std::vector<unsigned int> & nextseqtab) -> void {
    // refactoring: can data structures be marked as const?
    /* dump swarms */
    if (parameters.opt_mothur) {
      write_swarms_mothur_format(parameters, hashtable, nextseqtab);
    }
    else {
      write_swarms_default_format(parameters, hashtable, nextseqtab);
    }

    /* dump seeds in fasta format with sum of abundances */
    if (not parameters.opt_seeds.empty()) {
      write_representative_sequences(parameters, hashtable);
    }

    /* output swarm in uclust format */
    if (not parameters.opt_uclust_file.empty()) {
      write_swarms_uclust_format(parameters, hashtable, nextseqtab);
    }

    /* output internal structure to file */
    if (not parameters.opt_internal_structure.empty()) {
      write_structure_file(parameters, hashtable, nextseqtab);
    }

    /* output statistics to file */
    if (not parameters.opt_statistics_file.empty()) {
      write_stats_file(parameters, hashtable);
    }
  }
} // namespace


auto dereplicate(struct Parameters const & parameters) -> void
{
  const uint64_t dbsequencecount = db_getsequencecount();
  const uint64_t hashtablesize {compute_hashtable_size(dbsequencecount)};

  std::vector<struct bucket> hashtable(hashtablesize);
  /* alloc and init table of links to other sequences in cluster */
  std::vector<unsigned int> nextseqtab(dbsequencecount, 0);

  // dereplicate input sequences
  auto const stats = dereplicating(parameters, hashtable, nextseqtab);

  // sort by decreasing abundance
  sort_seeds(parameters, hashtable);
  release_unused_memory(hashtable, stats.swarmcount);

  output_results(parameters, hashtable, nextseqtab);

  std::fprintf(parameters.logfile, "\n");
  std::fprintf(parameters.logfile, "Number of swarms:  %" PRIu64 "\n",
               static_cast<uint64_t>(stats.swarmcount));
  std::fprintf(parameters.logfile, "Largest swarm:     %u\n", stats.maxsize);
  std::fprintf(parameters.logfile, "Heaviest swarm:    %" PRIu64 "\n", stats.maxmass);
}
