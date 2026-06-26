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

#include "algod1_network.h"
#include "../swarm.h"
#include "../db.h"
#include "variants.h"
#include "algod1_internal.h"
#include "bloom.h"
#include "hashtable.h"
#include "make_unique.h"
#include "progress.h"
#include "threads.h"
#include "view.h"
#include <algorithm>  // std::copy()
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <iterator>  // std::next()
#include <mutex>  // std::unique_lock
#include <vector>


namespace {

  inline auto find_variant_matches(struct Parameters const & parameters,
                                   Data const & data,
                                   Hashtable const & hash_table,
                                   BloomFilter const & bloom_a,
                                   unsigned int seed,
                                   struct var_s const & var,
                                   std::vector<unsigned int>& hits_data,
                                   unsigned int & hits_count) -> void
  {
    if (not bloom_a.get(var.hash)) {
      return;
    }

    /* compute hash and corresponding hash table index */

    auto index = hash_table.getindex(var.hash);

    /* find matching buckets */

    while (hash_table.is_occupied(index))
      {
        if (hash_table.compare_value(index, var.hash))
          {
            const auto amp = hash_table.get_data(index);

            /* avoid self */
            if (seed != amp) {
              if ((parameters.opt_no_cluster_breaking) or
                  (data.abundance(seed) >= data.abundance(amp)))
                {
                  auto const seed_seq = data.sequence_view(seed);
                  auto const amp_seq = data.sequence_view(amp);

                  if (check_variant(seed_seq, var, amp_seq))
                    {
                      hits_data[hits_count] = amp;
                      ++hits_count;
                      break;
                    }
                }
            }
          }
        index = hash_table.getnextindex(index);
      }
  }


  auto check_variants(struct Parameters const & parameters,
                      Data const & data,
                      Hashtable const & hash_table,
                      BloomFilter const & bloom_a,
                      unsigned int seed,
                      std::vector<struct var_s> & variant_list,
                      std::vector<unsigned int>& hits_data) -> unsigned int
  {
    auto hits_count = 0U;

    auto const seed_seq = data.sequence_view(seed);
    const auto hash = data.sequence_hash(seed);
    const auto variant_count = generate_variants(data.zobrist(), seed_seq, hash, variant_list);

    // variant_list is pre-sized to an upper bound; only the first
    // variant_count entries are valid for this call.
    auto const variants = View<var_s>{variant_list.data(), variant_count};
    for (auto const & var : variants) {
      find_variant_matches(parameters, data, hash_table, bloom_a, seed, var, hits_data, hits_count);
    }

    return hits_count;
  }


  auto network_thread(struct Parameters const & parameters,
                      Data const & data,
                      std::vector<struct ampinfo_s> & ampinfo_v,
                      Hashtable const & hash_table,
                      BloomFilter const & bloom_a,
                      struct Network_state & state,
                      Progress & progress) -> void
  {
    std::size_t const n_items = compute_microvariant_buffer_size(data.longest_sequence());

    std::vector<unsigned int> hits_data(n_items);
    std::vector<struct var_s> variant_list(n_items);

    const auto amplicons = data.sequence_count();
    std::unique_lock<std::mutex> lock(state.mutex);
    while (state.amp < amplicons)
      {
        const auto amp = state.amp;
        ++state.amp;
        progress.update(amp);

        lock.unlock();

        const auto hits_count = check_variants(parameters, data, hash_table, bloom_a, amp, variant_list, hits_data);
        lock.lock();

        auto & target_amplicon = ampinfo_v[amp];
        target_amplicon.link_start = state.count;
        target_amplicon.link_count = hits_count;

        while (state.count + hits_count > state.network_v.size()) {
          state.network_v.resize(state.network_v.size() + one_megabyte);
        }

        std::copy(hits_data.cbegin(),
                  std::next(hits_data.cbegin(), hits_count),
                  std::next(state.network_v.begin(),
                            static_cast<std::ptrdiff_t>(state.count)));
        state.count += hits_count;
      }
  }

} // namespace


auto build_amplicon_network(struct Parameters const & parameters,
                            Data const & data,
                            std::vector<struct ampinfo_s> & ampinfo_v,
                            struct Network_state & network_state) -> void
{
  /* d=1 hashtable and Bloom filter live in this function's scope so
     their backing storage is released before run_fastidious_pass()
     allocates its own fresh pair. */

  /* populate the d=1 hash table and Bloom filter with the amplicon
     hashes precomputed in db.cc */
  const auto amplicons = data.sequence_count();
  Hashtable hash_table;
  const auto hashtablesize = hash_table.allocate(amplicons);
  BloomFilter bloom_a(hashtablesize, amplicon_pattern_shift,
                      amplicon_n_hash_functions);

  Progress progress_hash("Building hashtable:", amplicons, parameters);

  for (auto k = 0U; k < amplicons; ++k)
    {
      hash_insert(data, hash_table, bloom_a, k);
      progress_hash.increment();
    }

  progress_hash.done();


  Progress progress_network("Building network: ", amplicons, parameters);
  {
    auto const network_tr = utils::make_unique<ThreadRunner>(
        parameters.opt_threads,
        [&parameters, &data, &ampinfo_v, &hash_table, &bloom_a, &network_state, &progress_network](uint64_t /*nth_thread*/) -> void {
          network_thread(parameters, data, ampinfo_v, hash_table, bloom_a, network_state, progress_network);
        });
    network_tr->run();
  }

  progress_network.done();
}
