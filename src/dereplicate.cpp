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

#include "dereplicate.hpp"
#include "swarm.hpp"
#include "db.hpp"
#include "utils/hashtable_size.hpp"
#include "utils/progress.hpp"
#include <algorithm>  // sort
#include <cassert>  // assert
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>
#include <cstdio>  // fputc()
#include <cstdlib>  // qsort()
#include <iterator>  // std::next
#include <vector>


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
    Progress const progress("Sorting:          ", 1, parameters);

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
      return lhs.seqno_first < rhs.seqno_first;
    };

    std::sort(hashtable.begin(), hashtable.end(), compare_seeds);
    progress.done();
  }


  auto release_unused_memory(std::vector<struct bucket>& hashtable,
                             int64_t const swarmcount) -> void {
    hashtable.erase(hashtable.begin() + swarmcount, hashtable.end());
    hashtable.shrink_to_fit();
    assert(hashtable.size() == static_cast<uint64_t>(swarmcount));
  }


  auto write_stats_file(struct Parameters const & parameters,
                        Data const & data,
                        std::vector<struct bucket> const & hashtable) -> void {
    Progress progress("Writing stats:    ", hashtable.size(), parameters);
    for (auto const & cluster: hashtable) {
      std::fprintf(parameters.statsfile.get(), "%u\t%" PRIu64 "\t", cluster.size, cluster.mass);
      data.fprint_id_noabundance(parameters.statsfile.get(), cluster.seqno_first, parameters.opt_usearch_abundance);
      std::fprintf(parameters.statsfile.get(), "\t%" PRIu64 "\t%u\t%u\t%u\n",
                   data.abundance(cluster.seqno_first),
                   cluster.singletons, 0U, 0U);
      progress.increment();
    }
    progress.done();
  }


  auto write_structure_file(struct Parameters const & parameters,
                            Data const & data,
                            std::vector<struct bucket> const & hashtable,
                            std::vector<unsigned int> const & nextseqtab) -> void {
    Progress progress("Writing structure:", hashtable.size(), parameters);
    auto counter = 0UL;

    for (auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;
      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          data.fprint_id_noabundance(parameters.internal_structure_file.get(), seed, parameters.opt_usearch_abundance);
          std::fprintf(parameters.internal_structure_file.get(), "\t");
          data.fprint_id_noabundance(parameters.internal_structure_file.get(), next_identical, parameters.opt_usearch_abundance);
          std::fprintf(parameters.internal_structure_file.get(), "\t%d\t%lu\t%d\n", 0, counter + 1, 0);
          next_identical = nextseqtab[next_identical];
        }
      ++counter;
      progress.update(counter);
    }
    progress.done();
  }


  auto write_swarms_uclust_format(struct Parameters const & parameters,
                                  Data const & data,
                                  std::vector<struct bucket> const & hashtable,
                                  std::vector<unsigned int> const & nextseqtab) -> void {
    Progress progress("Writing UCLUST:   ", hashtable.size(), parameters);
    auto counter = 0U;

    for (auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;

      std::fprintf(parameters.uclustfile.get(), "C\t%u\t%u\t*\t*\t*\t*\t*\t",
                   counter,
                   cluster.size);
      data.fprint_id(parameters.uclustfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile.get(), "\t*\n");

      std::fprintf(parameters.uclustfile.get(), "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                   counter,
                   data.sequence_view(seed).length);
      data.fprint_id(parameters.uclustfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(parameters.uclustfile.get(), "\t*\n");

      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          std::fprintf(parameters.uclustfile.get(),
                       "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                       counter,
                       data.sequence_view(next_identical).length,
                       100.0,
                       "=");
          data.fprint_id(parameters.uclustfile.get(), next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile.get(), "\t");
          data.fprint_id(parameters.uclustfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(parameters.uclustfile.get(), "\n");
          next_identical = nextseqtab[next_identical];
        }
      ++counter;
      progress.update(counter);
    }
    progress.done();
  }


  auto write_representative_sequences(struct Parameters const & parameters,
                                      Data const & data,
                                      std::vector<struct bucket> const & hashtable) -> void {
    Progress progress("Writing seeds:    ", hashtable.size(), parameters);
    for (auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;
      std::fprintf(parameters.seeds_file.get(), ">");
      data.fprint_id_with_new_abundance(parameters.seeds_file.get(), seed, cluster.mass, parameters.opt_usearch_abundance);
      std::fprintf(parameters.seeds_file.get(), "\n");
      data.fprintseq(parameters.seeds_file.get(), seed);
      progress.increment();
    }
    progress.done();
  }


  auto write_swarms_mothur_format(struct Parameters const & parameters,
                                  Data const & data,
                                  std::vector<struct bucket> const & hashtable,
                                  std::vector<unsigned int> const & nextseqtab) -> void {
    Progress progress("Writing swarms:   ", hashtable.size(), parameters);

#ifdef _WIN32
    std::fprintf(parameters.outfile.get(), "swarm_%" PRIu64 "\t%llu", parameters.opt_differences, hashtable.size());
#else
    std::fprintf(parameters.outfile.get(), "swarm_%" PRIu64 "\t%lu", parameters.opt_differences, hashtable.size());
#endif

    for (auto const & cluster: hashtable) {
      // print cluster seed
      auto const seed = cluster.seqno_first;
      std::fputc('\t', parameters.outfile.get());
      data.fprint_id(parameters.outfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      // print other cluster members
      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          std::fputc(',', parameters.outfile.get());
          data.fprint_id(parameters.outfile.get(), next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          next_identical = nextseqtab[next_identical];
        }

      progress.increment();
    }
    std::fputc('\n', parameters.outfile.get());

    progress.done();
  }


  auto write_swarms_default_format(struct Parameters const & parameters,
                                   Data const & data,
                                   std::vector<struct bucket> const & hashtable,
                                   std::vector<unsigned int> const & nextseqtab) -> void {
    static constexpr char sepchar {' '};
    Progress progress("Writing swarms:   ", hashtable.size(), parameters);

    for (auto const & cluster: hashtable) {
      // print cluster seed
      auto const seed = cluster.seqno_first;
      data.fprint_id(parameters.outfile.get(), seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      // print other cluster members
      auto next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          std::fputc(sepchar, parameters.outfile.get());
          data.fprint_id(parameters.outfile.get(), next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          next_identical = nextseqtab[next_identical];
        }
      std::fputc('\n', parameters.outfile.get());
      progress.increment();
    }

    progress.done();
  }


  auto dereplicating(struct Parameters const & parameters,
                     Data const & data,
                     std::vector<struct bucket> & hashtable,
                     std::vector<unsigned int> & nextseqtab)
    -> struct Stats
       {
         Progress progress("Dereplicating:    ", nextseqtab.size(), parameters);

         struct Stats stats;
         const uint64_t derep_hash_mask = hashtable.size() - 1;
         auto const & zobrist = data.zobrist();

         for (auto seqno = 0U; seqno < nextseqtab.size(); ++seqno)
           {
             auto const seq = data.sequence_view(seqno);

             /*
               Find free bucket or bucket for identical sequence.
               Make sure sequences are exactly identical
               in case of any hash collision.
               With 64-bit hashes, there is about 50% chance of a
               collision when the number of sequences is about 5e9.
             */

             auto const hash = zobrist.hash(seq);

             auto nth_bucket = hash & derep_hash_mask;
             auto * clusterp = &hashtable[nth_bucket];

             while ((clusterp->mass != 0U) and
                    ((clusterp->hash != hash) or
                     (seq.length != data.sequence_view(clusterp->seqno_first).length) or
                     not std::equal(seq.encoded.cbegin(), seq.encoded.cend(),
                                    data.sequence_view(clusterp->seqno_first).encoded.cbegin())
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

             auto const abundance = data.abundance(seqno);

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

             progress.update(seqno);
           }
         progress.done();

         return stats;
    }


  auto output_results(struct Parameters const & parameters,
                      Data const & data,
                      std::vector<struct bucket> const & hashtable,
                      std::vector<unsigned int> const & nextseqtab) -> void {
    /* dump swarms */
    if (parameters.opt_mothur) {
      write_swarms_mothur_format(parameters, data, hashtable, nextseqtab);
    }
    else {
      write_swarms_default_format(parameters, data, hashtable, nextseqtab);
    }

    /* dump seeds in fasta format with sum of abundances */
    if (not parameters.opt_seeds.empty()) {
      write_representative_sequences(parameters, data, hashtable);
    }

    /* output swarm in uclust format */
    if (not parameters.opt_uclust_file.empty()) {
      write_swarms_uclust_format(parameters, data, hashtable, nextseqtab);
    }

    /* output internal structure to file */
    if (not parameters.opt_internal_structure.empty()) {
      write_structure_file(parameters, data, hashtable, nextseqtab);
    }

    /* output statistics to file */
    if (not parameters.opt_statistics_file.empty()) {
      write_stats_file(parameters, data, hashtable);
    }
  }
} // namespace


auto dereplicate(struct Parameters const & parameters,
                 Data const & data) -> void
{
  const uint64_t dbsequencecount = data.sequence_count();
  const uint64_t hashtablesize {compute_hashtable_size(dbsequencecount)};

  std::vector<struct bucket> hashtable(hashtablesize);
  /* alloc and init table of links to other sequences in cluster */
  std::vector<unsigned int> nextseqtab(dbsequencecount, 0);

  // dereplicate input sequences
  auto const stats = dereplicating(parameters, data, hashtable, nextseqtab);

  // sort by decreasing abundance
  sort_seeds(parameters, hashtable);
  release_unused_memory(hashtable, stats.swarmcount);

  output_results(parameters, data, hashtable, nextseqtab);

  std::fprintf(parameters.logfile, "\n");
  std::fprintf(parameters.logfile, "Number of swarms:  %" PRIu64 "\n",
               static_cast<uint64_t>(stats.swarmcount));
  std::fprintf(parameters.logfile, "Largest swarm:     %u\n", stats.maxsize);
  std::fprintf(parameters.logfile, "Heaviest swarm:    %" PRIu64 "\n", stats.maxmass);
}
