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

#include "cpu_features.h"  // Cpu_features
#include <cstdint>  // uint64_t


// XOR + popcount over two qgramvectorbytes-long buffers.
//
// The implementation lives under src/arch/<isa>/qgram_compare.cc; the
// Makefile picks the right one via ARCH_DIR. Builds for an unsupported
// architecture fail to link (no static fallback — see #error in
// qgram.cc kept until a portable fallback is added).
//
// cpu_features is consulted only by the x86_64 path, which dispatches
// between SSE2 and POPCNT at run time. Other architectures ignore it.
auto compareqgramvectors(unsigned char const * lhs, unsigned char const * rhs,
                         Cpu_features const & cpu_features) -> uint64_t;

#endif
