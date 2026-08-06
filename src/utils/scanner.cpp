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

#include "../db.hpp"
#include "scanner.hpp"
#include "../search8.hpp"
#include "../search16.hpp"
#include "../swarm.hpp"
#include "cpu_features.hpp"  // Cpu_features
#include "memory_budget.hpp"  // require_ram
#include "score_matrix.hpp"
#include "search_data.hpp"  // Search_data, BYTE, WORD
#include "span.hpp"  // Span<uint64_t>
#include "threads.hpp"  // ThreadRunner
#include "view.hpp"  // View<uint64_t>
#include <array>
#include <cassert>  // assert()
#include <cstdint>  // int64_t, uint64_t
#include <mutex>  // std::lock_guard
#include <vector>

#ifndef NDEBUG
// C++17 refactoring: [[maybe_unused]]
#include <limits>
constexpr auto ullong_max = std::numeric_limits<unsigned long long int>::max();
#endif


namespace {

auto allocate_per_thread_search_data(std::vector<struct Search_data>& search_data_v,
                                     const uint64_t longestdbsequence) -> void {
  static constexpr auto one_kilobyte = 1024UL;
  static constexpr auto nt_per_uint64 = 32U;
  const uint64_t dirbuffersize = longestdbsequence * ((longestdbsequence + 3) / 4) * 4;

  // dir_array_v dominates and grows as O(L^2) in the longest sequence
  // length; fail early with a clear message instead of aborting inside
  // operator new (see memory_budget.hpp).
  static constexpr auto bytes_per_uint64 = uint64_t{8};
  const uint64_t per_thread_bytes =
      (dirbuffersize * bytes_per_uint64)                 // dir_array_v (dominant)
    + (longestdbsequence * nt_per_uint64)                // hearray_v
    + (longestdbsequence * 2 * sizeof(void *));          // qtable_v + qtable_w_v
  require_ram(per_thread_bytes,
             search_data_v.size(),
             "the pairwise-alignment buffers");

  for (auto & thread_data: search_data_v) {
    thread_data.qtable_v.resize(longestdbsequence);
    thread_data.qtable_w_v.resize(longestdbsequence);
    thread_data.dprofile_v.resize(2 * one_kilobyte);  // 4 * 16 * 32
    thread_data.dprofile_w_v.resize(1 * one_kilobyte);  // 4 * 2 * 8 * 32
    thread_data.hearray_v.resize(longestdbsequence * nt_per_uint64);
    thread_data.dir_array_v.resize(dirbuffersize);
  }
}


auto adjust_thread_number(const Bit_mode n_bits,
                          const uint64_t remaining_sequences,
                          uint64_t n_threads) -> uint64_t {
  static constexpr auto channels_8 = 8U;
  static constexpr auto channels_16 = 16U;
  const auto channels = (n_bits == Bit_mode::bits_16) ? channels_8 : channels_16;

  assert(remaining_sequences != 0);
  assert(n_threads != 0);
  assert((n_threads - 1) <= (ullong_max / channels_8));

  while (remaining_sequences <= (n_threads - 1) * channels) {
    --n_threads;
  }

  return n_threads;
}

// arguments: bits, length, thr
// static_assert(adjust_thread_number(Bit_mode::bits_8, 32, 10) == 2);
// static_assert(adjust_thread_number(Bit_mode::bits_8, 32,  3) == 2);
// static_assert(adjust_thread_number(Bit_mode::bits_8, 31,  2) == 2);
// static_assert(adjust_thread_number(Bit_mode::bits_8, 17,  2) == 2);
// static_assert(adjust_thread_number(Bit_mode::bits_8, 16,  2) == 1);
// static_assert(adjust_thread_number(Bit_mode::bits_8,  1,  2) == 1);
// static_assert(adjust_thread_number(Bit_mode::bits_8, 32,  1) == 1);
// static_assert(adjust_thread_number(Bit_mode::bits_16, 17, 10) == 3);
// static_assert(adjust_thread_number(Bit_mode::bits_16, 17,  3) == 3);
// static_assert(adjust_thread_number(Bit_mode::bits_16, 16,  3) == 2);
// static_assert(adjust_thread_number(Bit_mode::bits_16, 15,  2) == 2);
// static_assert(adjust_thread_number(Bit_mode::bits_16,  1,  3) == 1);
// static_assert(adjust_thread_number(Bit_mode::bits_16, 17,  1) == 1);

}  // namespace


Scanner::Scanner(struct Parameters const & parameters,
                 Data const & data)
  : data_(data),
    gapopen_(parameters.penalty_gapopen),
    gapextend_(parameters.penalty_gapextend),
    score_matrix_8_(create_score_matrix<unsigned char>(parameters.penalty_mismatch)),
    score_matrix_16_(create_score_matrix<unsigned short>(parameters.penalty_mismatch)),
    n_threads_(parameters.opt_threads.count()),
    search_data_v_(parameters.opt_threads.count()),
    threads_(parameters.opt_threads.count(),
             [this](uint64_t thread_id) -> void { worker_core(thread_id); }) {
  allocate_per_thread_search_data(search_data_v_, data.longest_sequence());

  Cpu_features const features {parameters.ssse3_present != 0,
                               parameters.sse41_present != 0,
                               parameters.popcnt_present != 0,};
  for (auto & thread_data : search_data_v_) {
    thread_data.cpu_features = features;
  }
}


auto Scanner::init(struct Search_data & thread_data) const -> void {
  static constexpr auto byte_multiplier = 64U;
  static constexpr auto word_multiplier = 32U;

  for (auto i = 0U; i < query_.length; ++i) {
    const auto nt_value = nucleotide_at(query_, i) + 1U;  // 1,  2,   3, or   4
    const auto byte_offset = byte_multiplier * nt_value;  // 1, 64, 128, or 192
    const auto word_offset = word_multiplier * nt_value;  // 1, 32,  64, or 128

    // refactoring: difficult to work directly on vectors (thread barrier)
    thread_data.qtable_v[i]   = &thread_data.dprofile_v[byte_offset];
    thread_data.qtable_w_v[i] = &thread_data.dprofile_w_v[word_offset];
  }
}


auto Scanner::chunk(struct Search_data & thread_data, const Bit_mode bits) -> void {
  assert(thread_data.target_count != 0);

  // The window this thread was handed by next_window(), computed once here
  // rather than as three unchecked pointer bumps: the subviews assert
  // their own bounds against the caller's arrays in debug builds.
  auto const first = thread_data.target_index;
  auto const count = thread_data.target_count;

  if (bits == Bit_mode::bits_16) {
    assert(gapopen_ <= std::numeric_limits<WORD>::max());
    assert(gapextend_ <= std::numeric_limits<WORD>::max());
    search16(data_.get(),
             thread_data,
             static_cast<WORD>(gapopen_),
             static_cast<WORD>(gapextend_),
             score_matrix_16_.data(),
             targets_.subview(first, count),
             scores_.subspan(first, count),
             diffs_.subspan(first, count),
             query_);
  } else {
    assert(gapopen_ <= std::numeric_limits<BYTE>::max());
    assert(gapextend_ <= std::numeric_limits<BYTE>::max());
    search8(data_.get(),
            thread_data,
            static_cast<BYTE>(gapopen_),
            static_cast<BYTE>(gapextend_),
            score_matrix_8_.data(),
            targets_.subview(first, count),
            scores_.subspan(first, count),
            diffs_.subspan(first, count),
            query_);
  }
}


auto Scanner::next_window() -> Scanner::Work_window {
  std::lock_guard<std::mutex> const lock(scan_mutex_);

  auto const listlength = targets_.size();
  if (next_ >= listlength) {
    return Work_window{};  // exhausted
  }

  const uint64_t chunksize =
    ((listlength - next_ + remainingchunks_ - 1) / remainingchunks_);
  Work_window const window {next_, chunksize};

  next_ += chunksize;
  --remainingchunks_;

  return window;
}


auto Scanner::worker_core(const uint64_t thread_id) -> void {
  auto & thread_data = search_data_v_[thread_id];
  init(thread_data);
  for (auto window = next_window(); not window.empty(); window = next_window()) {
    thread_data.target_index = window.first;
    thread_data.target_count = window.count;
    chunk(thread_data, bits_);
  }
}


auto Scanner::run(const uint64_t query_no,
                  View<uint64_t> const targets,
                  Span<uint64_t> const scores,
                  Span<uint64_t> const diffs,
                  const Bit_mode bits) -> void {
  assert(scores.size() == targets.size());
  assert(diffs.size() == targets.size());

  query_ = data_.get().sequence_view(query_no);

  next_ = 0;
  targets_ = targets;
  scores_ = scores;
  diffs_ = diffs;
  bits_ = bits;

  const auto thr = adjust_thread_number(bits, targets_.size(), n_threads_);

  remainingchunks_ = thr;

  if (thr == 1) {
    worker_core(0);
  }
  else {
    threads_.run();
  }
}
