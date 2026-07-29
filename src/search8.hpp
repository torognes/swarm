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

#include "utils/search_data.hpp"  // Search_data, BYTE
#include "utils/span.hpp"  // Span<uint64_t>
#include "utils/view.hpp"  // View<uint64_t>
#include <cstdint>


class Data;       // defined in db.hpp
struct Sequence;  // defined in db.hpp


// seqnos is the window of candidate amplicon indexes to align the query
// against; scores, diffs and alignmentlengths are the caller's output
// windows over the same candidates, so all four carry the same length.
auto search8(Data const & data,
             Search_data & search_data,
             BYTE gap_open_penalty,
             BYTE gap_extend_penalty,
             BYTE const * score_matrix,
             View<uint64_t> seqnos,
             Span<uint64_t> scores,
             Span<uint64_t> diffs,
             Span<uint64_t> alignmentlengths,
             Sequence const & query) -> void;
