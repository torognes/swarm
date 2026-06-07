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

#include "../db.h"  // Data (stored as reference_wrapper member)
#include "queryinfo.h"
#include "score_matrix.h"  // create_score_matrix, n_cells
#include "search_data.h"  // Search_data, BYTE, WORD
#include "threads.h"  // ThreadRunner
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <functional>  // std::reference_wrapper
#include <mutex>
#include <vector>


struct Parameters;  // defined in swarm.h


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

  // searches the query against listlength targets, writing scores,
  // diffs and alignment lengths back to the caller-owned arrays
  auto run(uint64_t query_no,
           uint64_t listlength,
           uint64_t * targets,
           uint64_t * scores,
           uint64_t * diffs,
           uint64_t * alignlengths,
           Bit_mode bits) -> void;

  // entry point for each worker thread (also called directly when a
  // single thread suffices)
  auto worker_core(uint64_t thread_id) -> void;

private:
  static constexpr std::size_t score_matrix_alignment {16};

  auto init(struct Search_data & thread_data) -> void;
  auto chunk(struct Search_data & thread_data, Bit_mode bits) -> void;
  auto getwork(uint64_t & countref, uint64_t & firstref) -> bool;

  std::reference_wrapper<Data const> data_;
  int64_t gapopen_ {0};
  int64_t gapextend_ {0};
  alignas(score_matrix_alignment)
    std::array<unsigned char, n_cells * n_cells> score_matrix_8_;
  alignas(score_matrix_alignment)
    std::array<unsigned short, n_cells * n_cells> score_matrix_16_;
  uint64_t n_threads_ {0};

  std::mutex scan_mutex_;
  struct queryinfo query_ {0, 0, nullptr};
  uint64_t master_next_ {0};
  uint64_t master_length_ {0};
  uint64_t remainingchunks_ {0};
  uint64_t * master_targets_ {nullptr};
  uint64_t * master_scores_ {nullptr};
  uint64_t * master_diffs_ {nullptr};
  uint64_t * master_alignlengths_ {nullptr};
  Bit_mode master_bits_ {Bit_mode::bits_16};

  std::vector<struct Search_data> search_data_v_;
  ThreadRunner threads_;  // last: its lambda touches the members above
};

#endif
