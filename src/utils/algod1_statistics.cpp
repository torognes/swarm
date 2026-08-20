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
#include "ceil_divide.hpp"  // ceil_divide
#include "system_memory.hpp"
#include "fatal.hpp"
#include "print_view.hpp"  // fprint, fprint_integer
#include "progress.hpp"
#include <algorithm>  // std::max()
#include <cassert>  // assert()
#include <cstdint>  // uint64_t
#include <cstdio>  // fprintf()
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
                            struct Bloom_demand const & demand) -> Bloom_geometry
{
  auto const nucleotides_in_small_clusters = demand.nucleotides;
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

  auto const memlimit = system_get_memlimit();
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
      uint64_t const memrest
        = (one_megabyte * parameters.opt_ceiling) - memused;
      auto const new_bits = n_bits_in_a_byte * memrest / (microvariants * nucleotides_in_small_clusters);
      if (new_bits < bits)
        {
          if (new_bits < 2) {
            fatal("Insufficient memory remaining for Bloom filter.");
          }
          fprint(parameters.logfile, "Reducing memory used for Bloom filter due to --ceiling option.\n");
          bits = new_bits;
          n_hash_functions = hash_functions_for(static_cast<unsigned int>(bits));
          bloom_length_in_bits = bloom_bits_for(bits);
        }
    }
  else if (nucleotides_in_small_clusters != 0)
    {
      // No --ceiling, so nothing bounded this filter: its size was a fixed
      // number of bits per nucleotide in a light cluster, and the only check
      // was the advisory warning below, after which the allocation was
      // attempted anyway. Not setting --ceiling was therefore what allowed
      // the filter to exceed memory and send the run into swap.
      //
      // The budget is what is free now, less what this phase still has to
      // allocate while the filter is alive (demand.headroom_bytes, counted by
      // its caller). memlimit rather than memtotal, so that a cgroup-capped
      // run is bounded by its cgroup and not by the host.
      //
      // Unlike the --ceiling branch this never fails: the user asked for no
      // limit, so a limit discovered here must not turn a run that works
      // today into one that refuses to start. When even the 2-bit floor does
      // not fit, the floor is used and the shortfall reported.
      auto const committed = memused + demand.headroom_bytes;
      auto const memrest = (memlimit > committed) ? (memlimit - committed) : uint64_t{0};
      auto const new_bits = n_bits_in_a_byte * memrest / (microvariants * nucleotides_in_small_clusters);
      if (new_bits < bits)
        {
          static constexpr uint64_t min_bits_per_entry {2};
          fprint(parameters.logfile,
                 "Reducing memory used for Bloom filter to fit available memory.\n");
          if (new_bits < min_bits_per_entry) {
            fprint(parameters.logfile,
                   "WARNING: available memory allows fewer than 2 bits per entry; "
                   "using 2. The fastidious pass will be slow.\n");
          }
          bits = std::max(new_bits, min_bits_per_entry);
          n_hash_functions = hash_functions_for(static_cast<unsigned int>(bits));
          bloom_length_in_bits = bloom_bits_for(bits);
        }
    }

  static constexpr uint64_t min_bloom_length_in_bits {64};  // at least 64 bits
  bloom_length_in_bits = std::max(bloom_length_in_bits, min_bloom_length_in_bits);

  // bloom_length is in bits (divide by 8 to get bytes)
  // bloom_length is guaranteed to be at least 64 (see code above)
  assert(bloom_length_in_bits != 0);  // safeguard for future changes
  assert(bloom_length_in_bits >= 64);
  // ceil_divide's precondition, named at the call site and not only inside
  // the template: the sum it forms must not overflow. bloom_length_in_bits
  // is a nucleotide count times 7 times at most 64 bits, so it cannot come
  // near this, but meeting the contract is the caller's job.
  assert(bloom_length_in_bits
         <= std::numeric_limits<uint64_t>::max() - (n_bits_in_a_byte - 1));

  // Rounded up, and computed once. The warning below and the returned
  // geometry are the same number of bytes, and used to convert it
  // separately, one rounding down and the other up -- while issue 174 was
  // about a bits-to-bytes conversion that rounded the wrong way.
  auto const n_bytes = ceil_divide<uint64_t>(bloom_length_in_bits, n_bits_in_a_byte);

  if (memused + n_bytes > memlimit)
    {
      fprint(parameters.logfile, "WARNING: Memory usage will probably exceed total amount of memory available.\n");
      fprint(parameters.logfile, "Try to reduce memory footprint using the --bloom-bits or --ceiling options.\n");
    }

  fprint(parameters.logfile, "Bloom filter: bits=");
  fprint_integer(parameters.logfile, bits);
  fprint(parameters.logfile, ", m=");
  fprint_integer(parameters.logfile, bloom_length_in_bits);
  fprint(parameters.logfile, ", k=");
  fprint_integer(parameters.logfile, n_hash_functions);
  // the size in MB is the one field that needs a double formatted, so it
  // keeps its fprintf (see the same note in the uclust writers)
  static_cast<void>(std::fprintf(parameters.logfile, ", size=%.1fMB\n",
                                 static_cast<double>(bloom_length_in_bits) / (n_bits_in_a_byte * one_megabyte)));


  Bloom_geometry geom;
  geom.n_bytes = n_bytes;
  geom.n_hash_functions = n_hash_functions;
  return geom;
}


auto log_swarm_summary(struct Parameters const & parameters,
                       Overall_stats const & overall_stats) -> void
{
  fprint(parameters.logfile, '\n');

  fprint(parameters.logfile, "Number of swarms:  ");
  fprint_integer(parameters.logfile, overall_stats.swarmcount_adjusted);
  fprint(parameters.logfile, '\n');

  fprint(parameters.logfile, "Largest swarm:     ");
  fprint_integer(parameters.logfile, overall_stats.largest);
  fprint(parameters.logfile, '\n');

  fprint(parameters.logfile, "Max generations:   ");
  fprint_integer(parameters.logfile, overall_stats.maxgen);
  fprint(parameters.logfile, '\n');
}
