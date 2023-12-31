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
#include "util.h"
#include "utils/hashtable_size.h"
#include "utils/nt_codec.h"
#include "zobrist.h"
#include <cassert>  // assert
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstdint>
#include <cstdio>  // fputc()
#include <cstdlib>  // qsort()
#include <cstring>  // memcmp
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


struct bucket
{
  uint64_t hash = 0;
  unsigned int seqno_first = 0;
  unsigned int seqno_last = 0;
  uint64_t mass = 0;
  unsigned int size = 0;
  unsigned int singletons = 0;
};


// refactoring: std::stable_sort uses input order (no need to use seqno)
auto derep_compare(const void * a, const void * b) -> int
{
  const auto * x = static_cast<const struct bucket *>(a);
  const auto * y = static_cast<const struct bucket *>(b);
  int status {0};

  /* highest abundance first, otherwise keep order */

  if (x->mass < y->mass) {
    status = +1;
  }
  else if (x->mass > y->mass) {
    status = -1;
  }
  else {
    if (x->seqno_first < y->seqno_first) {
      status = -1;
    }
    else if (x->seqno_first > y->seqno_first) {
      status = +1;
    }
    else {
      status = 0;
    }
  }

  return status;
}


auto write_stats_file(struct Parameters const & parameters,
                      std::vector<struct bucket> & hashtable) -> void {
  progress_init("Writing stats:    ", hashtable.size());
  for(auto i = 0ULL; i < hashtable.size(); i++)
    {
      struct bucket * sp = hashtable.data() + i;
      std::fprintf(statsfile, "%u\t%" PRIu64 "\t", sp->size, sp->mass);
      fprint_id_noabundance(statsfile, sp->seqno_first, parameters.opt_usearch_abundance);
      std::fprintf(statsfile, "\t%" PRIu64 "\t%u\t%u\t%u\n",
              db_getabundance(sp->seqno_first),
              sp->singletons, 0U, 0U);
      progress_update(i);
    }
  progress_done();
}


auto write_structure_file(const uint64_t swarmcount,
                          struct Parameters const & parameters,
                          struct bucket * hashtable,
                          unsigned int * nextseqtab) -> void {
  progress_init("Writing structure:", swarmcount);

  for(uint64_t i = 0; i < swarmcount; i++)
    {
      struct bucket * sp = hashtable + i;
      const uint64_t seed = sp->seqno_first;
      unsigned int next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          fprint_id_noabundance(internal_structure_file, seed, parameters.opt_usearch_abundance);
          std::fprintf(internal_structure_file, "\t");
          fprint_id_noabundance(internal_structure_file, next_identical, parameters.opt_usearch_abundance);
          std::fprintf(internal_structure_file, "\t%d\t%" PRIu64 "\t%d\n", 0, i + 1, 0);
          next_identical = nextseqtab[next_identical];
        }
      progress_update(i);
    }
  progress_done();
}


auto write_swarms_uclust_format(const uint64_t swarmcount,
                                struct Parameters const & parameters,
                                struct bucket * hashtable,
                                unsigned int * nextseqtab) -> void {
  progress_init("Writing UCLUST:   ", swarmcount);

  for(auto swarmid = 0U; swarmid < swarmcount; swarmid++)
    {
      struct bucket * bp = hashtable + swarmid;

      const unsigned int seed = bp->seqno_first;

      std::fprintf(uclustfile, "C\t%u\t%u\t*\t*\t*\t*\t*\t",
              swarmid,
              bp->size);
      fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(uclustfile, "\t*\n");

      std::fprintf(uclustfile, "S\t%u\t%u\t*\t*\t*\t*\t*\t",
              swarmid,
              db_getsequencelen(seed));
      fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
      std::fprintf(uclustfile, "\t*\n");

      unsigned int next_identical = nextseqtab[seed];
      while (next_identical != 0U)
        {
          std::fprintf(uclustfile,
                  "H\t%u\t%u\t%.1f\t+\t0\t0\t%s\t",
                  swarmid,
                  db_getsequencelen(next_identical),
                  100.0,
                  "=");
          fprint_id(uclustfile, next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(uclustfile, "\t");
          fprint_id(uclustfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
          std::fprintf(uclustfile, "\n");
          next_identical = nextseqtab[next_identical];
        }

      progress_update(swarmid + 1);
    }
  progress_done();
}


auto write_representative_sequences(struct Parameters const & parameters,
                                    std::vector<struct bucket> const & hashtable) -> void {
  progress_init("Writing seeds:    ", hashtable.size());
  auto counter = 0U;
  for(auto const & cluster: hashtable) {
    const unsigned int seed = cluster.seqno_first;
    std::fprintf(fp_seeds, ">");
    fprint_id_with_new_abundance(fp_seeds, seed, cluster.mass, parameters.opt_usearch_abundance);
    std::fprintf(fp_seeds, "\n");
    db_fprintseq(fp_seeds, seed);
    ++counter;
    progress_update(counter);
  }
  progress_done();
}


auto write_swarms_mothur_format(struct Parameters const & parameters,
                                std::vector<struct bucket> const & hashtable,
                                std::vector<unsigned int> const & nextseqtab) -> void {
  progress_init("Writing swarms:   ", hashtable.size());
  std::fprintf(outfile, "swarm_%" PRId64 "\t%" PRIu64, parameters.opt_differences, hashtable.size());
  auto counter = 0U;

  for(auto const & cluster: hashtable) {
    // print cluster seed
    const unsigned int seed = cluster.seqno_first;
    std::fputc('\t', outfile);
    fprint_id(outfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

    // print other cluster members
    unsigned int next_identical = nextseqtab[seed];
    while (next_identical != 0U)
      {
        std::fputc(',', outfile);
        fprint_id(outfile, next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        next_identical = nextseqtab[next_identical];
      }

    ++counter;
    progress_update(counter);
  }
    std::fputc('\n', outfile);

  progress_done();
}


auto write_swarms_default_format(struct Parameters const & parameters,
                                 std::vector<struct bucket> const & hashtable,
                                 std::vector<unsigned int> const & nextseqtab) -> void {
  progress_init("Writing swarms:   ", hashtable.size());
  auto counter = 0U;

  for(auto const & cluster: hashtable) {
    // print cluster seed
    const unsigned int seed = cluster.seqno_first;
    fprint_id(outfile, seed, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

    // print other cluster members
    unsigned int next_identical = nextseqtab[seed];
    while (next_identical != 0U)
      {
        std::fputc(sepchar, outfile);
        fprint_id(outfile, next_identical, parameters.opt_usearch_abundance, parameters.opt_append_abundance);
        next_identical = nextseqtab[next_identical];
      }
    std::fputc('\n', outfile);
    ++counter;
    progress_update(counter);
  }

  progress_done();
}


auto dereplicate(struct Parameters const & parameters) -> void
{
  const uint64_t dbsequencecount = db_getsequencecount();
  const uint64_t hashtablesize {compute_hashtable_size(dbsequencecount)};
  const uint64_t derep_hash_mask = hashtablesize - 1;

  std::vector<struct bucket> hashtable(hashtablesize);

  uint64_t swarmcount = 0;
  int64_t tmp_swarmcount = 0;
  uint64_t maxmass = 0;
  unsigned int maxsize = 0;

  /* alloc and init table of links to other sequences in cluster */
  std::vector<unsigned int> nextseqtab(dbsequencecount, 0);

  progress_init("Dereplicating:    ", dbsequencecount);

  for(auto i = 0U; i < dbsequencecount; i++)
    {
      const unsigned int seqlen = db_getsequencelen(i);
      char * seq = db_getsequence(i);

      /*
        Find free bucket or bucket for identical sequence.
        Make sure sequences are exactly identical
        in case of any hash collision.
        With 64-bit hashes, there is about 50% chance of a
        collision when the number of sequences is about 5e9.
      */

      const uint64_t hash = zobrist_hash(reinterpret_cast<unsigned char *>(seq),
                                         seqlen);

      uint64_t j = hash & derep_hash_mask;
      struct bucket * bp = hashtable.data() + j;

      while (((bp->mass) != 0U) and
             ((bp->hash != hash) or
              (seqlen != db_getsequencelen(bp->seqno_first)) or
              (std::memcmp(seq,
                      db_getsequence(bp->seqno_first),
                      nt_bytelength(seqlen)) != 0)))
        {
          ++bp;
          ++j;
          if (bp >= hashtable.data() + hashtablesize) // wrap around the table if we reach the end
            {
              bp = hashtable.data();
              j = 0;
            }
        }

      const uint64_t abundance = db_getabundance(i);

      if ((bp->mass) != 0U)
        {
          /* at least one identical sequence already */
          nextseqtab[bp->seqno_last] = i;
        }
      else
        {
          /* no identical sequences yet, start a new cluster */
          ++swarmcount;
          ++tmp_swarmcount;
          bp->hash = hash;
          bp->seqno_first = i;
          bp->size = 0;
          bp->singletons = 0;
        }

      ++bp->size;
      bp->seqno_last = i;
      bp->mass += abundance;

      if (abundance == 1) {
        ++bp->singletons;
      }

      if (bp->mass > maxmass) {
        maxmass = bp->mass;
      }

      if (bp->size > maxsize) {
        maxsize = bp->size;
      }

      progress_update(i);
    }
  progress_done();

  progress_init("Sorting:          ", 1);
  std::qsort(hashtable.data(), hashtablesize, sizeof(bucket), derep_compare);
  hashtable.erase(hashtable.begin() + tmp_swarmcount, hashtable.end());
  hashtable.shrink_to_fit();  // release unused memory
  assert(hashtable.size() == swarmcount);
  progress_done();


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
  if (uclustfile != nullptr) {
    write_swarms_uclust_format(swarmcount, parameters, hashtable.data(), nextseqtab.data());
  }

  /* output internal structure to file */
  if (not parameters.opt_internal_structure.empty()) {
    write_structure_file(swarmcount, parameters, hashtable.data(), nextseqtab.data());
  }

  /* output statistics to file */
  if (statsfile != nullptr) {
    write_stats_file(parameters, hashtable);
  }

  std::fprintf(logfile, "\n");
  std::fprintf(logfile, "Number of swarms:  %" PRIu64 "\n", swarmcount);
  std::fprintf(logfile, "Largest swarm:     %u\n", maxsize);
  std::fprintf(logfile, "Heaviest swarm:    %" PRIu64 "\n", maxmass);
}
