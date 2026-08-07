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

#ifndef SWARM_UTILS_QGRAM_COMPARE_H
#define SWARM_UTILS_QGRAM_COMPARE_H

#include "cpu_features.hpp"  // Cpu_features
#include "qgram_array.hpp"  // Qgram_vector
#include <cstdint>  // uint64_t


// XOR + popcount over two qgram vectors.
//
// The implementation lives under src/arch/<isa>/qgram_compare.cpp; the
// Makefile picks the right one via ARCH_DIR. Builds for an unsupported
// architecture fail to link (no static fallback — see #error in
// qgram.cpp kept until a portable fallback is added).
//
// Qgram_vector rather than unsigned char const *: every SIMD kernel below
// reads its argument with aligned loads, and Qgram_vector is the type that
// carries the alignas(16) making that legal. Taking the bytes instead
// erased it at the seam, leaving the requirement stated in a comment that
// no assert could check -- a pointer's alignment is not recoverable from
// unsigned char const *.
//
// cpu_features is consulted only by the x86_64 path, which dispatches
// between SSE2 and POPCNT at run time. Other architectures ignore it.
auto compareqgramvectors(Qgram_vector const & lhs, Qgram_vector const & rhs,
                         Cpu_features const & cpu_features) -> uint64_t;

#endif
