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

#ifndef SWARM_UTILS_ALGOD1_FASTIDIOUS_H
#define SWARM_UTILS_ALGOD1_FASTIDIOUS_H


/* Fastidious phase for the d=1 algorithm: enumerate the second-
   generation microvariants of heavy-swarm seeds, look for matches
   against the light-swarm amplicons via a dedicated Bloom filter, and
   graft the matching light swarms onto their heavy parents. */

#include "algod1_internal.hpp"
#include <vector>

struct Parameters;  // defined in swarm.hpp


auto run_fastidious_pass(struct Parameters const & parameters,
                         Data const & data,
                         unsigned int swarmcount,
                         std::vector<struct ampinfo_s> & ampinfo_v,
                         std::vector<struct swarminfo_s> & swarminfo_v,
                         Overall_stats & overall_stats) -> void;

#endif  // SWARM_UTILS_ALGOD1_FASTIDIOUS_H
