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

#include "qgram.hpp"
#include "../db.hpp"
#include "../swarm.hpp"
#include "ceil_divide.hpp"
#include "cpu_features.hpp"
#include "progress.hpp"
#include "qgram_array.hpp"
#include "qgram_compare.hpp"  // compareqgramvectors (per-arch impl under arch/<isa>/)
#include "span.hpp"  // Span<uint64_t>
#include "threads.hpp"
#include "view.hpp"  // View<uint64_t>
#include <algorithm>  // std::transform
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <limits>
#include <vector>


namespace {

  auto findqgrams(Sequence const & sequence,
                  Qgram_vector & qgramvector) noexcept -> void
  {
    /* set qgram bit vector by xoring occurrences of qgrams in sequence */

    static constexpr unsigned int max_range {7};

    qgramvector.fill(0);

    auto const seqlen = sequence.length;
    uint64_t qgram {0};
    unsigned int position {0};

    while ((position < qgramlength - 1) and (position < seqlen)) {
      qgram = (qgram << 2U) | nucleotide_at(sequence, position);
      ++position;
    }

    while (position < seqlen) {
      qgram = (qgram << 2U) | nucleotide_at(sequence, position);
      assert((qgram & max_range) <= 7);
      auto const index = (qgram >> 3U) & (qgramvectorbytes - 1);
      qgramvector[index] ^= static_cast<unsigned char>(1U << (qgram & max_range));
      ++position;
    }
  }


  auto build_qgram_store(struct Parameters const & parameters,
                         Data const & data) -> Qgram_store
  {
    auto const n_sequences = data.sequence_count();
    Qgram_store store(n_sequences);

    Progress progress_qg("Find qgram vects: ", n_sequences, parameters);
    // increment() rather than update(counter): counter is a 0-based index,
    // so handing it over left the bar one milestone short of n_sequences
    for (auto counter = 0U; counter < n_sequences; ++counter) {
      findqgrams(data.sequence_view(counter), store[counter]);
      progress_qg.increment();
    }
    progress_qg.done();
    return store;
  }


  inline auto qgram_diff(Qgram_store const & store,
                         uint64_t const seqno_a, uint64_t const seqno_b,
                         Cpu_features const & cpu_features) noexcept -> uint64_t
  {
    assert(seqno_a < store.size());
    assert(seqno_b < store.size());
    uint64_t const diffqgrams = compareqgramvectors(store[seqno_a],
                                                    store[seqno_b],
                                                    cpu_features);
    // Each mismatch flips up to 2*qgramlength bits in the qgram XOR
    // vector (q bits leave, q bits arrive). Dividing the bit-difference
    // count by that upper bound gives a lower bound on edit distance.
    static constexpr uint64_t max_bits_per_mismatch {qgramlength + qgramlength};
    return ceil_divide(diffqgrams, max_bits_per_mismatch);  // mindiff
  }

}  // namespace


QgramDiffer::QgramDiffer(struct Parameters const & parameters,
                         Data const & data)
  : store_(build_qgram_store(parameters, data)),
    cpu_features_{
      parameters.ssse3_present != 0,
      parameters.sse41_present != 0,
      parameters.popcnt_present != 0,
    },
    thread_info_v_(parameters.opt_threads.count()),
    threads_(parameters.opt_threads.count(),
             [this](uint64_t nth_thread) -> void {
               worker(nth_thread);
             })
{ }


auto QgramDiffer::worker(uint64_t const nth_thread) const noexcept -> void
{
  assert(nth_thread < thread_info_v_.size());
  auto const & tip = thread_info_v_[nth_thread];

  auto const seed = tip.seed;
  auto const difflist = tip.difflist;

  // one input shape or the other, never both, see thread_info_s. Not
  // "exactly one": the last cluster of a run finds the pool already
  // exhausted, so both are empty and the transform below is a no-op.
  assert(tip.amplist.empty() or tip.poollist.empty());

  // an empty chunk takes this branch too, and transforms nothing
  if (tip.poollist.empty()) {
    // one distance per candidate, so the chunk's two halves agree in length
    assert(difflist.size() == tip.amplist.size());
    std::transform(tip.amplist.cbegin(), tip.amplist.cend(), difflist.begin(),
                   [this, seed](uint64_t const candidate) noexcept -> uint64_t {
                     return qgram_diff(store_, seed, candidate, cpu_features_);
                   });
    return;
  }

  assert(difflist.size() == tip.poollist.size());
  std::transform(tip.poollist.cbegin(), tip.poollist.cend(), difflist.begin(),
                 [this, seed](struct ampliconinfo_s const & candidate) noexcept -> uint64_t {
                   return qgram_diff(store_, seed, candidate.ampliconid, cpu_features_);
                 });
}


template <typename Assign>
auto QgramDiffer::distribute_and_run(uint64_t const seed,
                                     uint64_t const listlen,
                                     Assign assign) -> void
{
  static constexpr auto single_threaded_threshold = std::numeric_limits<uint8_t>::max();
  if (listlen <= single_threaded_threshold)
    {
      auto & tip = thread_info_v_[0];
      tip.seed = seed;
      assign(tip, std::size_t{0}, static_cast<std::size_t>(listlen));
      worker(0);
      return;
    }

  std::size_t offset {0};
  auto listrest = listlen;
  auto thrrest = thread_info_v_.size();

  /* distribute work */
  for (auto & tip: thread_info_v_) {
      auto const chunk = ceil_divide(listrest, thrrest);

      tip.seed = seed;
      assign(tip, offset, static_cast<std::size_t>(chunk));

      offset += chunk;
      listrest -= chunk;
      --thrrest;
    }

  threads_.run();
}


auto QgramDiffer::fast(uint64_t const seed,
                       View<uint64_t> const amplist,
                       Span<uint64_t> const difflist) -> void
{
  assert(difflist.size() == amplist.size());
  distribute_and_run(seed, amplist.size(),
                     [amplist, difflist](thread_info_s & tip,
                                         std::size_t const offset,
                                         std::size_t const chunk) -> void {
                       tip.amplist = amplist.subview(offset, chunk);
                       tip.poollist = View<struct ampliconinfo_s>{};
                       tip.difflist = difflist.subspan(offset, chunk);
                     });
}


auto QgramDiffer::fast_over_pool(uint64_t const seed,
                                 View<struct ampliconinfo_s> const poollist,
                                 Span<uint64_t> const difflist) -> void
{
  assert(difflist.size() == poollist.size());
  distribute_and_run(seed, poollist.size(),
                     [poollist, difflist](thread_info_s & tip,
                                          std::size_t const offset,
                                          std::size_t const chunk) -> void {
                       tip.amplist = View<uint64_t>{};
                       tip.poollist = poollist.subview(offset, chunk);
                       tip.difflist = difflist.subspan(offset, chunk);
                     });
}
