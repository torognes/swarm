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

#ifndef SWARM_UTILS_SIMD_ALIGNMENT_H
#define SWARM_UTILS_SIMD_ALIGNMENT_H

#include <cstddef>  // std::size_t


// The vector registers swarm's kernels use are 128 bits wide, so 16 bytes is
// both how much one load moves and the boundary that load has to start on.
// Five buffers are accessed that way, all taking the figure from here rather
// than repeating it:
//
//   Qgram_vector  (qgram_array.hpp)  read by the qgram_compare kernels
//   He_block      (search_data.hpp)  the HE array, read by search8/search16
//   Dseq_8/16     (search_data.hpp)  read by dprofile_shuffle8/16
//   Dprofile_8/16 (search_data.hpp)  written by the score-profile builders
//   score_matrix_8_ and _16_ (scanner.hpp)  read by dprofile_fill8/16
//
// The first four say so in their own types; the score matrices say so at
// their declaration instead, because create_score_matrix() returns the
// plain std::array that Score_matrix_8/16 alias.
//
// Widening this would not be enough on its own to move to 256-bit vectors:
// the kernels name their vector type per architecture, and the score-matrix
// offsets in dprofile_fill8/16 assume 16-byte steps.
constexpr std::size_t simd_vector_bytes {16};

#endif  // SWARM_UTILS_SIMD_ALIGNMENT_H
