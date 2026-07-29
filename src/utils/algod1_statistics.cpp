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

#include "algod1_statistics.hpp"
#include "../swarm.hpp"
#include "algod1_internal.hpp"
#include "system_memory.hpp"
#include "fatal.hpp"
#include "progress.hpp"
#include <algorithm>  // std::max()
#include <cassert>  // assert()
#include <cinttypes>  // macro PRIu64
#include <cstdint>  // uint64_t
#include <cstdio>  // fprintf(), fputc(), fputs()
#include <vector>  // std::vector

#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
#include <limits>
constexpr auto uint_max = std::numeric_limits<unsigned int>::max();
#endif


auto count_cluster_stats(struct Parameters const & parameters,
                         unsigned int const amplicon_count,
                         std::vector<struct swarminfo_s> const & swarminfo_v) -> Cluster_stats
{
  Cluster_stats stats;

  Progress progress_count("Counting amplicons in heavy and light swarms",
                          swarminfo_v.size(), parameters);

  for (auto const & swarm_info : swarminfo_v)
    {
      if (swarm_info.mass < parameters.opt_boundary)
        {
          stats.amplicons_in_small_clusters += swarm_info.size;
          stats.nucleotides_in_small_clusters += swarm_info.sumlen;
          ++stats.small_clusters;
        }
      progress_count.increment();
    }
  progress_count.done();

  stats.amplicons_in_large_clusters = amplicon_count - stats.amplicons_in_small_clusters;
  stats.large_clusters = swarminfo_v.size() - stats.small_clusters;

  return stats;
}


auto compute_bloom_geometry(struct Parameters const & parameters,
                            uint64_t const nucleotides_in_small_clusters) -> Bloom_geometry
{
  /* m: total size of Bloom filter in bits */
  /* k: number of hash functions (n_hash_functions) */
  /* n: number of entries in the bloom filter */
  /* here: k=11 and m/n=18, that is 16 bits/entry */

  static constexpr auto microvariants = 7U;
  static constexpr auto n_bits_in_a_byte = 8U;
  static constexpr double hash_functions_per_bit {4.0 / 10};
  static constexpr double natural_log_of_2 {0.693147181};  // C++26 refactoring: std::log(2.0)
  static_assert(hash_functions_per_bit <= natural_log_of_2, "upper limit is log(2)");
  assert(parameters.opt_bloom_bits <= uint_max);
  assert(parameters.opt_bloom_bits <= 64);  // larger than expected
  assert(parameters.opt_bloom_bits >= 2);  // smaller than expected

  // bits_value is unsigned int (not uint64_t) to avoid a risky
  // uint64 to double conversion warning in the multiplication
  auto const hash_functions_for = [](unsigned int const bits_value) -> unsigned int {
    return std::max(static_cast<unsigned int>(hash_functions_per_bit * bits_value), 1U);
  };
  auto const bloom_bits_for = [nucleotides_in_small_clusters](uint64_t const bits_value) -> uint64_t {
    return nucleotides_in_small_clusters * microvariants * bits_value;
  };

  auto bits = parameters.opt_bloom_bits;

  // int64_t n_hash_functions = int(bits * std::log(2.0));    /* 16 bits -> 11 hash functions */
  // auto n_hash_functions = unsigned int(hash_functions_per_bit * bits); /* 6 */
  auto n_hash_functions = hash_functions_for(static_cast<unsigned int>(bits));

  auto bloom_length_in_bits = bloom_bits_for(bits);

  auto const memtotal = system_get_memtotal();
  auto const memused = system_get_memused();

  // nucleotides_in_small_clusters guards the division below: it is the
  // denominator of new_bits. The caller only reaches this function when
  // there is at least one light swarm, so it is non-zero in practice;
  // the explicit check keeps the ceiling adjustment safe regardless.
  if ((parameters.opt_ceiling != 0) and (nucleotides_in_small_clusters != 0))
    {
      if (parameters.opt_ceiling * one_megabyte < memused)
        {
          fatal("Memory ceiling for Bloom filter is too low.");
        }
      assert(memused < one_megabyte * parameters.opt_ceiling);
      const uint64_t memrest
        = (one_megabyte * parameters.opt_ceiling) - memused;
      auto const new_bits = n_bits_in_a_byte * memrest / (microvariants * nucleotides_in_small_clusters);
      if (new_bits < bits)
        {
          if (new_bits < 2) {
            fatal("Insufficient memory remaining for Bloom filter.");
          }
          static_cast<void>(std::fputs("Reducing memory used for Bloom filter due to --ceiling option.\n", parameters.logfile));
          bits = new_bits;
          n_hash_functions = hash_functions_for(static_cast<unsigned int>(bits));
          bloom_length_in_bits = bloom_bits_for(bits);
        }
    }

  static constexpr uint64_t min_bloom_length_in_bits {64};  // at least 64 bits
  bloom_length_in_bits = std::max(bloom_length_in_bits, min_bloom_length_in_bits);

  if (memused + (bloom_length_in_bits / n_bits_in_a_byte) > memtotal)
    {
      static_cast<void>(std::fputs("WARNING: Memory usage will probably exceed total amount of memory available.\n", parameters.logfile));
      static_cast<void>(std::fputs("Try to reduce memory footprint using the --bloom-bits or --ceiling options.\n", parameters.logfile));
    }

  std::fprintf(parameters.logfile,
               "Bloom filter: bits=%" PRIu64 ", m=%" PRIu64 ", k=%u, size=%.1fMB\n",
               bits, bloom_length_in_bits, n_hash_functions, static_cast<double>(bloom_length_in_bits) / (n_bits_in_a_byte * one_megabyte));


  // bloom_length is in bits (divide by 8 to get bytes)
  // bloom_length is guaranteed to be at least 64 (see code above)
  assert(bloom_length_in_bits != 0);  // safeguard for future changes
  assert(bloom_length_in_bits >= 64);

  Bloom_geometry geom;
  geom.n_bytes = ((bloom_length_in_bits - 1) / n_bits_in_a_byte) + 1;
  geom.n_hash_functions = n_hash_functions;
  return geom;
}


auto log_swarm_summary(struct Parameters const & parameters,
                       Overall_stats const & overall_stats) -> void
{
  static_cast<void>(std::fputc('\n', parameters.logfile));
  std::fprintf(parameters.logfile, "Number of swarms:  %" PRIu64 "\n", overall_stats.swarmcount_adjusted);
  std::fprintf(parameters.logfile, "Largest swarm:     %u\n", overall_stats.largest);
  std::fprintf(parameters.logfile, "Max generations:   %u\n", overall_stats.maxgen);
}
