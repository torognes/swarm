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
#include "simd_alignment.hpp"  // simd_vector_bytes
#include <array>
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
