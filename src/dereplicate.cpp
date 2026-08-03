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
#include "utils/chain_range.hpp"  // Chain_range
#include "utils/hashtable_size.hpp"
#include "utils/print_view.hpp"  // fprint, fprint_integer
#include "utils/progress.hpp"
#include <algorithm>  // sort
#include <cassert>  // assert
#include <cinttypes>  // macro PRIu64
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fprintf(), fputc(), fputs()
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


  // The sequences identical to a cluster's seed are a chain through
  // nextseqtab, ending at 0 (see utils/chain_range.hpp).
  struct Next_identical
  {
    static auto next(std::vector<unsigned int> const & nextseqtab,
                     unsigned int const seqno) -> unsigned int {
      return nextseqtab[seqno];
    }
  };

  using Identical_chain = Chain_range<std::vector<unsigned int>, Next_identical>;

  // The copies of `seed`, seed excluded: every consumer reports the seed
  // itself first and then walks the copies.
  inline auto identical_copies_of(std::vector<unsigned int> const & nextseqtab,
                                  unsigned int const seed) -> Identical_chain {
    return Identical_chain{nextseqtab, nextseqtab[seed], 0U};
  }


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
    auto const cluster_count = static_cast<std::size_t>(swarmcount);
    // resize() lengthens a vector that is too short, so shrinking is the
    // caller's invariant rather than something resize() enforces:
    // swarmcount counts occupied buckets and cannot exceed the table.
    assert(cluster_count <= hashtable.size());
    hashtable.resize(cluster_count);
    hashtable.shrink_to_fit();
  }


  auto write_stats_file(struct Parameters const & parameters,
                        Data const & data,
                        std::vector<struct bucket> const & hashtable) -> void {
    Progress progress("Writing stats:    ", hashtable.size(), parameters);
    for (auto const & cluster: hashtable) {
      std::fprintf(parameters.statsfile.get(), "%u\t%" PRIu64 "\t", cluster.size, cluster.mass);
      fprint_id_noabundance(parameters.statsfile.get(), data.info(cluster.seqno_first), parameters.opt_usearch_abundance);
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
    auto counter = uint64_t{0};

    auto * const structure_file = parameters.internal_structure_file.get();

    for (auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;
      for (auto const next_identical : identical_copies_of(nextseqtab, seed))
        {
          fprint_id_noabundance(structure_file, data.info(seed), parameters.opt_usearch_abundance);
          fprint(structure_file, '\t');
          fprint_id_noabundance(structure_file, data.info(next_identical), parameters.opt_usearch_abundance);
          // Columns 3 and 5 are the number of differences and the number of
          // steps from the seed. At d = 0 the two amplicons are identical, so
          // both are zero for every pair and go in as literal text.
          fprint(structure_file, "\t0\t");
          fprint_integer(structure_file, counter + 1);
          fprint(structure_file, "\t0\n");
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
      fprint_id(parameters.uclustfile.get(), data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      static_cast<void>(std::fputs("\t*\n", parameters.uclustfile.get()));

      std::fprintf(parameters.uclustfile.get(), "S\t%u\t%u\t*\t*\t*\t*\t*\t",
                   counter,
                   data.sequence_view(seed).length);
      fprint_id(parameters.uclustfile.get(), data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      static_cast<void>(std::fputs("\t*\n", parameters.uclustfile.get()));

      for (auto const next_identical : identical_copies_of(nextseqtab, seed))
        {
          std::fprintf(parameters.uclustfile.get(),
                       "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                       counter,
                       data.sequence_view(next_identical).length,
                       100.0,
                       "=");
          fprint_id(parameters.uclustfile.get(), data.info(next_identical), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          static_cast<void>(std::fputc('\t', parameters.uclustfile.get()));
          fprint_id(parameters.uclustfile.get(), data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          static_cast<void>(std::fputc('\n', parameters.uclustfile.get()));
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
    auto * const seeds_file = parameters.seeds_file.get();
    // one scratch buffer for the whole file
    Sequence_printer const sequence_printer {data.longest_sequence()};

    for (auto const & cluster: hashtable) {
      auto const seed = cluster.seqno_first;
      static_cast<void>(std::fputc('>', seeds_file));
      fprint_id_with_new_abundance(seeds_file, data.info(seed), cluster.mass,
                                   parameters.opt_usearch_abundance);
      static_cast<void>(std::fputc('\n', seeds_file));
      sequence_printer.print(seeds_file, data.sequence_view(seed));
      progress.increment();
    }
    progress.done();
  }


  auto write_swarms_mothur_format(struct Parameters const & parameters,
                                  Data const & data,
                                  std::vector<struct bucket> const & hashtable,
                                  std::vector<unsigned int> const & nextseqtab) -> void {
    Progress progress("Writing swarms:   ", hashtable.size(), parameters);

    uint64_t const number_of_clusters {hashtable.size()};
    std::fprintf(parameters.outfile.get(), "swarm_%" PRIu64 "\t%" PRIu64,
                 parameters.opt_differences, number_of_clusters);

    for (auto const & cluster: hashtable) {
      // print cluster seed
      auto const seed = cluster.seqno_first;
      static_cast<void>(std::fputc('\t', parameters.outfile.get()));
      fprint_id(parameters.outfile.get(), data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      // print other cluster members
      for (auto const next_identical : identical_copies_of(nextseqtab, seed))
        {
          static_cast<void>(std::fputc(',', parameters.outfile.get()));
          fprint_id(parameters.outfile.get(), data.info(next_identical), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        }

      progress.increment();
    }
    static_cast<void>(std::fputc('\n', parameters.outfile.get()));

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
      fprint_id(parameters.outfile.get(), data.info(seed), parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      // print other cluster members
      for (auto const next_identical : identical_copies_of(nextseqtab, seed))
        {
          static_cast<void>(std::fputc(sepchar, parameters.outfile.get()));
          fprint_id(parameters.outfile.get(), data.info(next_identical), parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        }
      static_cast<void>(std::fputc('\n', parameters.outfile.get()));
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

  static_cast<void>(std::fputc('\n', parameters.logfile));
  std::fprintf(parameters.logfile, "Number of swarms:  %" PRIu64 "\n",
               static_cast<uint64_t>(stats.swarmcount));
  std::fprintf(parameters.logfile, "Largest swarm:     %u\n", stats.maxsize);
  std::fprintf(parameters.logfile, "Heaviest swarm:    %" PRIu64 "\n", stats.maxmass);
}
