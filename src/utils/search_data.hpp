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

#ifndef SWARM_UTILS_SEARCH_DATA_H
#define SWARM_UTILS_SEARCH_DATA_H

#include "cpu_features.hpp"
#include "score_matrix.hpp"  // n_cells
#include "simd_alignment.hpp"  // simd_vector_bytes
#include <array>
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <vector>


using BYTE = unsigned char;
using WORD = unsigned short;

// alignas: both kernels hand hearray_v to their inner loop as a
// VECTORTYPE * and read it with aligned loads (__m128i on x86_64,
// uint8x16_t on aarch64, vector unsigned char on ppc), so the buffer has to
// start on a 16-byte boundary. A std::vector<BYTE> satisfied that only by
// accident: operator new returns storage aligned for any fundamental type,
// which is exactly 16 on these targets, so the requirement was met with no
// margin and nothing in the code said it existed. An over-aligned element
// type states it, and lets the kernels cast from a type whose own alignment
// already covers the target instead of from bare bytes. Deriving from
// std::array keeps data() and size() available.
struct alignas(simd_vector_bytes) He_block : std::array<BYTE, simd_vector_bytes> {};

// The two read-only inputs of the score-profile builders, named. Their
// destination is dprofile_v below, which already had a name.
//
// dprofile_fill8/16, dprofile_shuffle8/16 and dispatch_dprofile8/16 used to
// take all three as same-family pointers in a row -- dprofile, score_matrix,
// dseq -- where nothing but the argument order said which was which, and no
// extent was visible to the callee.

// The substitution scores, as returned by create_score_matrix() and held by
// Scanner, which over-aligns them: dprofile_fill8/16 loads 16 bytes at a
// time from an offset that is a multiple of 32 (see utils/scanner.hpp).
using Score_matrix_8  = std::array<BYTE, n_cells * n_cells>;
using Score_matrix_16 = std::array<WORD, n_cells * n_cells>;

// SIMD lanes at each search width: 16 at 8 bits, 8 at 16 bits. search8.cpp
// and search16.cpp assert their own 'channels' against these.
//
// Spelt out rather than suffixed with the width, because the count and the
// width are each other's suffix here: 8 bits gives 16 channels. Everything
// else in this family -- Score_matrix_8, Dseq_8, search8, dprofile_fill8 --
// is suffixed with the width.
constexpr std::size_t channels_at_8_bits {simd_vector_bytes / sizeof(BYTE)};
constexpr std::size_t channels_at_16_bits {simd_vector_bytes / sizeof(WORD)};

// The staging buffer holding the next block of database nucleotides: one
// SIMD vector's worth of bytes per channel, laid out as cdepth blocks of
// 'channels' bytes each (see utils/dseq_fill.hpp).
using Dseq_8  = std::array<BYTE, channels_at_8_bits * simd_vector_bytes>;
using Dseq_16 = std::array<BYTE, channels_at_16_bits * simd_vector_bytes>;

struct Search_data
{
  std::vector<BYTE *> qtable_v;
  std::vector<WORD *> qtable_w_v;

  std::vector<BYTE> dprofile_v;
  std::vector<WORD> dprofile_w_v;

  std::vector<He_block> hearray_v;  // sized in blocks, not bytes
  std::vector<uint64_t> dir_array_v;

  uint64_t target_count = 0;
  uint64_t target_index = 0;

  Cpu_features cpu_features {false, false, false};
};

#endif
