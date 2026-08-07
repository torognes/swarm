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

#ifndef SWARM_UTILS_QGRAM_ARRAY_H
#define SWARM_UTILS_QGRAM_ARRAY_H

#include "simd_alignment.hpp"  // simd_vector_bytes
#include <array>
#include <vector>


// 128 bytes = 1,024 bits, one bit per possible 5-mer (4^5 = 1,024).
// qgramvectorbytes is derived from qgramlength below.
constexpr unsigned int qgramlength     {5};
constexpr unsigned int qgramvectorbytes {(1U << (2 * qgramlength)) / 8};

// alignas(16): the per-architecture qgram_compare kernels read these with
// aligned SIMD loads (_mm_load_si128 on x86_64, and the NEON/VMX
// equivalents), so each element must start on a 16-byte boundary. A bare
// std::array has alignment 1; deriving an over-aligned type makes
// std::vector<Qgram_vector> place its buffer -- and, since the 128-byte
// element size is a multiple of 16, every element -- on a 16-byte
// boundary. Mirrors the alignas(16) on the score matrices in scanner.hpp.
struct alignas(simd_vector_bytes) Qgram_vector : std::array<unsigned char, qgramvectorbytes> {};
using Qgram_store  = std::vector<Qgram_vector>;

#endif
