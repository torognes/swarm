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

#ifndef SWARM_SEARCH16_H
#define SWARM_SEARCH16_H


#include "utils/search_data.hpp"  // Search_data, WORD
#include "utils/span.hpp"  // Span<uint64_t>
#include "utils/view.hpp"  // View<uint64_t>
#include <cstdint>


class Data;       // defined in db.hpp
struct Sequence;  // defined in db.hpp


// seqnos is the window of candidate amplicon indexes to align the query
// against; scores and diffs are the caller's output windows over the same
// candidates, so all three carry the same length.
auto search16(Data const & data,
              Search_data & search_data,
              WORD gap_open_penalty,
              WORD gap_extend_penalty,
              WORD const * score_matrix,
              View<uint64_t> seqnos,
              Span<uint64_t> scores,
              Span<uint64_t> diffs,
              Sequence const & query) -> void;

#endif  // SWARM_SEARCH16_H
