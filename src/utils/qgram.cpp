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
#include "nt_codec.hpp"
#include "threads.hpp"
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <iterator>  // std::next
#include <limits>
#include <vector>


namespace {

  auto findqgrams(Sequence const & sequence,
                  Qgram_vector & qgramvector) noexcept -> void
  {
    /* set qgram bit vector by xoring occurrences of qgrams in sequence */

    static constexpr unsigned int max_range {7};

    qgramvector.fill(0);

    auto const * const seq = sequence.encoded.data();
    auto const seqlen = sequence.length;
    uint64_t qgram {0};
    unsigned int position {0};

    while ((position < qgramlength - 1) and (position < seqlen)) {
      qgram = (qgram << 2U) | nt_extract(seq, position);
      ++position;
    }

    while (position < seqlen) {
      qgram = (qgram << 2U) | nt_extract(seq, position);
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
    for (auto counter = 0U; counter < n_sequences; ++counter) {
      findqgrams(data.sequence_view(counter), store[counter]);
      progress_qg.update(counter);
    }
    progress_qg.done();
    return store;
  }


  inline auto qgram_diff(Qgram_store const & store,
                         uint64_t seqno_a, uint64_t seqno_b,
                         Cpu_features const & cpu_features) noexcept -> uint64_t
  {
    assert(seqno_a < store.size());
    assert(seqno_b < store.size());
    const uint64_t diffqgrams = compareqgramvectors(store[seqno_a].data(),
                                                    store[seqno_b].data(),
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
    thread_info_v_(parameters.opt_threads),
    threads_(parameters.opt_threads,
             [this](uint64_t nth_thread) -> void {
               worker(nth_thread);
             })
{ }


auto QgramDiffer::worker(uint64_t const nth_thread) const noexcept -> void
{
  assert(nth_thread < thread_info_v_.size());
  auto const & tip = thread_info_v_[nth_thread];

  const auto seed = tip.seed;
  const auto listlen = tip.listlen;
  assert(listlen <= std::numeric_limits<std::ptrdiff_t>::max());
  const auto listlen_signed = static_cast<int64_t>(listlen);
  auto const * amplist = tip.amplist;
  auto * difflist = tip.difflist;

  for (auto i = 0LL; i < listlen_signed; ++i) {
    auto & target_diff = *std::next(difflist, i);
    auto const target_amplicon = *std::next(amplist, i);
    target_diff = qgram_diff(store_, seed, target_amplicon, cpu_features_);
  }
}


auto QgramDiffer::fast(uint64_t seed,
                       uint64_t const listlen,
                       std::vector<uint64_t> const & amplist,
                       std::vector<uint64_t> & difflist) -> void
{
  assert(listlen <= amplist.size());
  assert(listlen <= difflist.size());
  static constexpr auto single_threaded_threshold = std::numeric_limits<uint8_t>::max();
  if (listlen <= single_threaded_threshold)
    {
      auto & tip = thread_info_v_[0];
      tip.seed = seed;
      tip.listlen = listlen;
      tip.amplist = amplist.data();
      tip.difflist = difflist.data();
      worker(0);
    }
  else
    {
      auto const * next_amplist = amplist.data();
      auto * next_difflist = difflist.data();
      auto listrest = listlen;
      auto thrrest = thread_info_v_.size();

      /* distribute work */
      for (auto & tip: thread_info_v_) {
          auto const chunk = (listrest + thrrest - 1) / thrrest;
          assert(chunk <= std::numeric_limits<std::ptrdiff_t>::max());
          auto const chunk_signed = static_cast<int64_t>(chunk);

          tip.seed = seed;
          tip.listlen = chunk;
          tip.amplist = next_amplist;
          tip.difflist = next_difflist;

          next_amplist = std::next(next_amplist, chunk_signed);
          next_difflist = std::next(next_difflist, chunk_signed);
          listrest -= chunk;
          --thrrest;
        }

      threads_.run();
    }
}
