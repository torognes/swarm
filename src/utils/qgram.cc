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

#include "qgram.h"
#include "../db.h"
#include "../swarm.h"
#include "cpu_features.h"
#include "progress.h"
#include "qgram_array.h"
#include "qgram_compare.h"  // compareqgramvectors (per-arch impl under arch/<isa>/)
#include "nt_codec.h"
#include "threads.h"
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <cstring>  // memset
#include <iterator>  // std::next
#include <limits>
#include <vector>


namespace {
  // The q in q-gram, and the intermediate bit count derived from it.
  // qgramvectorbytes (= qgramvectorbits / 8) is the public constant in
  // utils/qgram_array.h; the static_assert below catches any drift.
  constexpr unsigned int qgramlength     {5};
  constexpr unsigned int qgramvectorbits {1U << (2 * qgramlength)};  // 1,024
  static_assert(qgramvectorbytes == qgramvectorbits / 8,
                "qgramvectorbytes must equal 4^qgramlength / 8");


  auto findqgrams(char const * seq, uint64_t seqlen,
                  unsigned char * qgramvector) -> void
  {
    /* set qgram bit vector by xoring occurrences of qgrams in sequence */

    static constexpr unsigned int max_range {7};

    std::memset(qgramvector, 0, qgramvectorbytes);

    uint64_t qgram {0};
    unsigned int position {0};

    while ((position < qgramlength - 1) and (position < seqlen)) {
      qgram = (qgram << 2U) | nt_extract(seq, position);
      ++position;
    }

    while (position < seqlen) {
      qgram = (qgram << 2U) | nt_extract(seq, position);
      assert((qgram & max_range) <= 7);
      assert(((qgram >> 3U) & (qgramvectorbytes - 1)) <= std::numeric_limits<std::ptrdiff_t>::max());
      auto const signed_position = static_cast<std::ptrdiff_t>((qgram >> 3U) & (qgramvectorbytes - 1));
      auto & target_qgram = *std::next(qgramvector, signed_position);
      target_qgram ^= static_cast<unsigned char>(1U << (qgram & max_range));
      ++position;
    }
  }
}  // namespace


namespace {

inline auto db_getqgramvector(Qgram_store const & store, uint64_t const seqno) -> unsigned char const *
{
  assert(seqno < store.size());
  return store[seqno].data();
}


auto build_qgram_store(struct Parameters const & parameters,
                       Data const & data) -> Qgram_store
{
  auto const n_sequences = data.sequence_count();
  Qgram_store store(n_sequences);

  Progress progress_qg("Find qgram vects: ", n_sequences, parameters);
  for (auto counter = 0U; counter < n_sequences; ++counter) {
    auto const seq = data.sequence_view(counter);
    findqgrams(seq.encoded.data(), seq.length, store[counter].data());
    progress_qg.update(counter);
  }
  progress_qg.done();
  return store;
}


inline auto qgram_diff(Qgram_store const & store,
                       uint64_t seqno_a, uint64_t seqno_b,
                       Cpu_features const & cpu_features) -> uint64_t
{
  const uint64_t diffqgrams = compareqgramvectors(db_getqgramvector(store, seqno_a),
                                                  db_getqgramvector(store, seqno_b),
                                                  cpu_features);
  return (diffqgrams + (2ULL * qgramlength) - 1) / (2ULL * qgramlength);  // mindiff
}

}  // namespace


QgramDiffer::QgramDiffer(struct Parameters const & parameters,
                         Data const & data)
  : store_(build_qgram_store(parameters, data)),
    cpu_features_{
      parameters.ssse3_present != 0,
      parameters.sse41_present != 0,
      parameters.popcnt_present != 0
    },
    thread_info_v_(static_cast<uint64_t>(parameters.opt_threads)),
    threads_(static_cast<std::size_t>(parameters.opt_threads),
             [this](uint64_t nth_thread) -> void {
               worker(nth_thread);
             })
{ }


auto QgramDiffer::worker(uint64_t const nth_thread) const -> void
{
  auto const & tip = *std::next(thread_info_v_.begin(), static_cast<std::ptrdiff_t>(nth_thread));

  const auto seed = tip.seed;
  const auto listlen = tip.listlen;
  assert(listlen <= std::numeric_limits<std::ptrdiff_t>::max());
  const auto listlen_signed = static_cast<int64_t>(listlen);
  auto * amplist = tip.amplist;
  auto * difflist = tip.difflist;

  for (auto i = 0LL; i < listlen_signed; ++i) {
    auto & target_diff = *std::next(difflist, i);
    auto const target_amplicon = *std::next(amplist, i);
    target_diff = qgram_diff(store_, seed, target_amplicon, cpu_features_);
  }
}


auto QgramDiffer::fast(uint64_t seed,
                       std::vector<uint64_t> const & amplist,
                       std::vector<uint64_t> & difflist) -> void
{
  assert(amplist.size() <= difflist.size());
  auto const listlen = amplist.size();
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
