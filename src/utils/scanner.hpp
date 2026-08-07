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

#ifndef SWARM_UTILS_SCANNER_H
#define SWARM_UTILS_SCANNER_H

#include "../db.hpp"  // Data (stored as reference_wrapper member), Sequence
#include "score_matrix.hpp"  // create_score_matrix, n_cells
#include "search_data.hpp"  // Search_data, BYTE, WORD
#include "simd_alignment.hpp"  // simd_vector_bytes
#include "span.hpp"  // Span<uint64_t>
#include "thread_count.hpp"  // ThreadCount
#include "threads.hpp"  // ThreadRunner
#include "view.hpp"  // View<uint64_t>
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <functional>  // std::reference_wrapper
#include <mutex>
#include <vector>


struct Parameters;  // defined in swarm.hpp


// Selects which SIMD search kernel processes a chunk: search8 packs the
// scores into 8-bit channels, search16 into 16-bit channels (chosen when
// 8 bits cannot hold the score, see set_bit_mode).
enum struct Bit_mode : std::uint8_t { bits_8, bits_16 };


// Aligns the score matrices used by the SIMD search kernels on 8- and
// 16-bit channels respectively, fanning the per-query work out over a
// pool of worker threads.
class Scanner {
public:
  Scanner(struct Parameters const & parameters,
          Data const & data);

  // Non-copyable, non-movable: the ThreadRunner's lambda captures
  // `this`, so the object must keep a stable address.
  Scanner(Scanner const &) = delete;
  Scanner(Scanner &&) = delete;
  auto operator=(Scanner const &) -> Scanner & = delete;
  auto operator=(Scanner &&) -> Scanner & = delete;
  ~Scanner() = default;

  // Searches the query against the amplicons listed in 'targets',
  // writing scores and diffs back to the caller-owned windows. All three
  // must be the same length: they are one candidate list and its two
  // result columns, indexed in lockstep.
  auto run(uint64_t query_no,
           View<uint64_t> targets,
           Span<uint64_t> scores,
           Span<uint64_t> diffs,
           Bit_mode bits) -> void;

  // entry point for each worker thread (also called directly when a
  // single thread suffices)
  auto worker_core(uint64_t thread_id) -> void;

private:
  // A thread's share of the target list: `count` entries starting at
  // `first`. Returned as one value rather than written through two
  // adjacent uint64_t references, which a caller could fill in either
  // order. An empty window (count == 0) means the list is exhausted.
  //
  // No default member initializers: under C++11 they would make this a
  // non-aggregate, and next_window() returns Work_window{first, count}.
  // Work_window{} value-initializes both members, which is the empty
  // window; there is no uninitialized declaration of this type.
  // C++14 refactoring: add {0} initializers, aggregates may have them
  //
  // Both members stay public, and misc-non-private-member-variables-in-
  // classes reports both, because empty() makes this a class with a member
  // function rather than a plain aggregate. Private members would need a
  // two-uint64_t constructor, which is the swappable pair of arguments the
  // paragraph above says this type exists to prevent.
  struct Work_window {
    uint64_t first;
    uint64_t count;

    constexpr auto empty() const noexcept -> bool { return count == 0; }
  };

  auto init(struct Search_data & thread_data) const -> void;
  auto chunk(struct Search_data & thread_data, Bit_mode bits) -> void;
  auto next_window() -> Work_window;

  std::reference_wrapper<Data const> data_;
  int64_t gapopen_ {0};
  int64_t gapextend_ {0};
  alignas(simd_vector_bytes)
    std::array<unsigned char, n_cells * n_cells> score_matrix_8_;
  alignas(simd_vector_bytes)
    std::array<unsigned short, n_cells * n_cells> score_matrix_16_;
  ThreadCount n_threads_ {};

  std::mutex scan_mutex_;
  Sequence query_ {};
  uint64_t next_ {0};
  uint64_t remainingchunks_ {0};
  // one candidate list and its two result columns; targets_.size() is
  // the list length that used to be tracked separately in length_.
  // No {} initializer: View and Span default-construct empty through
  // their own member initializers, so one here would be redundant.
  View<uint64_t> targets_;
  Span<uint64_t> scores_;
  Span<uint64_t> diffs_;
  Bit_mode bits_ {Bit_mode::bits_16};

  std::vector<struct Search_data> search_data_v_;
  ThreadRunner threads_;  // last: its lambda touches the members above
};

#endif
