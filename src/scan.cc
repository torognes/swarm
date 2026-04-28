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

#include "db.h"
#include "scan.h"
#include "search8.h"
#include "search16.h"
#include "swarm.h"
#include "utils/nt_codec.h"
#include "utils/search_data.h"
#include "utils/score_matrix.h"
#include <cassert>  // assert()
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <iterator>
#include <mutex>  // std::lock_guard
#include <vector>

#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
#include <limits>
constexpr auto ullong_max = std::numeric_limits<unsigned long long int>::max();
#endif


auto allocate_per_thread_search_data(std::vector<struct Search_data>& search_data_v,
                                     const uint64_t longestdbsequence) -> void
{
  static constexpr auto one_kilobyte = 1024UL;
  static constexpr auto nt_per_uint64 = 32U;
  const uint64_t dirbuffersize = longestdbsequence * ((longestdbsequence + 3) / 4) * 4;

  for(auto& thread_data: search_data_v) {
    thread_data.qtable_v.resize(longestdbsequence);
    thread_data.qtable_w_v.resize(longestdbsequence);
    thread_data.dprofile_v.resize(2 * one_kilobyte);  // 4 * 16 * 32
    thread_data.dprofile_w_v.resize(1 * one_kilobyte);  // 4 * 2 * 8 * 32
    thread_data.hearray_v.resize(longestdbsequence * nt_per_uint64);
    thread_data.dir_array_v.resize(dirbuffersize);
  }
}


auto search_init(struct Search_data & thread_data,
                 struct queryinfo const & query) -> void
{
  static constexpr auto byte_multiplier = 64U;
  static constexpr auto word_multiplier = 32U;

  for(auto i = 0U; i < query.len; ++i)
  {
    const auto nt_value = nt_extract(query.seq, i) + 1U;  // 1,  2,   3, or   4
    const auto byte_offset = byte_multiplier * nt_value;  // 1, 64, 128, or 192
    const auto word_offset = word_multiplier * nt_value;  // 1, 32,  64, or 128

    // refactoring: difficult to work directly on vectors (thread barrier)
    thread_data.qtable_v[i]   = &thread_data.dprofile_v[byte_offset];
    thread_data.qtable_w_v[i] = &thread_data.dprofile_w_v[word_offset];
  }
}


auto search_chunk(struct Parameters const & parameters,
                  Data const & data,
                  struct Search_data & thread_data,
                  struct Search_state const & state,
                  const int64_t bits) -> void
{
  static constexpr auto sixteen_bytes = 16;
  alignas(sixteen_bytes) static auto score_matrix_8 = create_score_matrix<unsigned char>(parameters.penalty_mismatch);
  alignas(sixteen_bytes) static auto score_matrix_16 = create_score_matrix<unsigned short>(parameters.penalty_mismatch);
  static constexpr auto bit_mode_16 = 16U;
  assert(thread_data.target_index <= std::numeric_limits<std::ptrdiff_t>::max());
  auto const target_index = static_cast<std::ptrdiff_t>(thread_data.target_index);

  assert(thread_data.target_count != 0);
  assert((bits == bit_mode_16) or (bits == bit_mode_16 / 2));

  if (bits == bit_mode_16) {
    assert(parameters.penalty_gapopen <= std::numeric_limits<WORD>::max());
    assert(parameters.penalty_gapextend <= std::numeric_limits<WORD>::max());
    search16(data,
             thread_data.qtable_w_v,
             static_cast<WORD>(parameters.penalty_gapopen),
             static_cast<WORD>(parameters.penalty_gapextend),
             score_matrix_16.data(),
             thread_data.dprofile_w_v,
             reinterpret_cast<WORD *>(thread_data.hearray_v.data()),
             thread_data.target_count,
             std::next(state.master_targets, target_index),
             std::next(state.master_scores, target_index),
             std::next(state.master_diffs, target_index),
             std::next(state.master_alignlengths, target_index),
             state.query.seq,
             static_cast<uint64_t>(state.query.len),
             thread_data.dir_array_v,
             thread_data.cpu_features);
  } else {
    assert(parameters.penalty_gapopen <= std::numeric_limits<BYTE>::max());
    assert(parameters.penalty_gapextend <= std::numeric_limits<BYTE>::max());
    search8(data,
            thread_data.qtable_v,
            static_cast<BYTE>(parameters.penalty_gapopen),
            static_cast<BYTE>(parameters.penalty_gapextend),
            score_matrix_8.data(),
            thread_data.dprofile_v,
            thread_data.hearray_v.data(),
            thread_data.target_count,
            std::next(state.master_targets, target_index),
            std::next(state.master_scores, target_index),
            std::next(state.master_diffs, target_index),
            std::next(state.master_alignlengths, target_index),
            state.query.seq,
            static_cast<uint64_t>(state.query.len),
            thread_data.dir_array_v,
            thread_data.cpu_features);
  }
}


auto search_getwork(struct Search_state & state,
                    uint64_t & countref, uint64_t & firstref) -> bool
{
  // countref = how many sequences to search
  // firstref = index into master_targets/scores/diffs where thread should start

  bool status {false};

  std::lock_guard<std::mutex> const lock(state.scan_mutex);

  if (state.master_next < state.master_length)
    {
      const uint64_t chunksize =
        ((state.master_length - state.master_next + state.remainingchunks - 1) / state.remainingchunks);

      countref = chunksize;
      firstref = state.master_next;

      state.master_next += chunksize;
      --state.remainingchunks;
      status = true;
    }

  return status;
}


auto search_worker_core(struct Parameters const & parameters,
                        Data const & data,
                        const int64_t thread_id, struct Search_state & state) -> void {
  auto & thread_data = *std::next(state.search_data, thread_id);
  search_init(thread_data, state.query);
  while(search_getwork(state, thread_data.target_count, thread_data.target_index)) {
    search_chunk(parameters, data, thread_data, state, state.master_bits);
  }
}


auto adjust_thread_number(const int n_bits,
                          const uint64_t remaining_sequences,
                          uint64_t n_threads) -> uint64_t {
  static constexpr auto channels_8 = 8U;
  static constexpr auto channels_16 = 16U;
  static constexpr auto bit_mode_16 = 16U;  // refactoring: should be an enum class
  const auto channels = (n_bits == bit_mode_16) ? channels_8 : channels_16;

  assert(remaining_sequences != 0);
  assert(n_threads != 0);
  assert((n_threads - 1) <= (ullong_max / channels_8));
  assert((n_bits == bit_mode_16) or (n_bits == bit_mode_16 / 2));

  while (remaining_sequences <= (n_threads - 1) * channels) {
    --n_threads;
  }

  return n_threads;
}

// arguments: bits, master_length, thr
// static_assert(adjust_thread_number( 8, 32, 10) == 2);
// static_assert(adjust_thread_number( 8, 32,  3) == 2);
// static_assert(adjust_thread_number( 8, 31,  2) == 2);
// static_assert(adjust_thread_number( 8, 17,  2) == 2);
// static_assert(adjust_thread_number( 8, 16,  2) == 1);
// static_assert(adjust_thread_number( 8,  1,  2) == 1);
// static_assert(adjust_thread_number( 8, 32,  1) == 1);
// static_assert(adjust_thread_number(16, 17, 10) == 3);
// static_assert(adjust_thread_number(16, 17,  3) == 3);
// static_assert(adjust_thread_number(16, 16,  3) == 2);
// static_assert(adjust_thread_number(16, 15,  2) == 2);
// static_assert(adjust_thread_number(16,  1,  3) == 1);
// static_assert(adjust_thread_number(16, 17,  1) == 1);


auto search_do(struct Parameters const & parameters,
               Data const & data,
               struct Search_state & state,
               const uint64_t query_no,
               const uint64_t listlength,
               uint64_t * targets,
               uint64_t * scores,
               uint64_t * diffs,
               uint64_t * alignlengths,
               const int bits,
               ThreadRunner * search_threads) -> void
{
  auto query_len = 0U;
  state.query.qno = query_no;
  auto const & info = data.info(query_no);
  state.query.seq = info.seq;
  query_len = info.seqlen;
  state.query.len = query_len;

  state.master_next = 0;
  state.master_length = listlength;
  state.master_targets = targets;
  state.master_scores = scores;
  state.master_diffs = diffs;
  state.master_alignlengths = alignlengths;
  state.master_bits = bits;

  const auto thr =
    adjust_thread_number(bits,
                         state.master_length,
                         static_cast<uint64_t>(parameters.opt_threads));

  state.remainingchunks = thr;

  if (thr == 1) {
    search_worker_core(parameters, data, 0, state);
  }
  else {
    search_threads->run();
  }
}


auto search_begin(struct Parameters const & parameters,
                  Data const & data,
                  struct Search_state & state,
                  std::vector<struct Search_data> & search_data_v) -> void
{
  state.search_data = search_data_v.data();

  allocate_per_thread_search_data(search_data_v, data.longest_sequence());

  for (auto & thread_data : search_data_v) {
    thread_data.cpu_features.ssse3 = (parameters.ssse3_present != 0);
    thread_data.cpu_features.sse41 = (parameters.sse41_present != 0);
    thread_data.cpu_features.popcnt = (parameters.popcnt_present != 0);
  }
}


auto search_end(struct Search_state & state) -> void
{
  /* finish and clean up worker threads */

  state.search_data = nullptr;
}
