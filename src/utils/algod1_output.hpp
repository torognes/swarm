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

#ifndef SWARM_UTILS_ALGOD1_OUTPUT_H
#define SWARM_UTILS_ALGOD1_OUTPUT_H


/* Result writers for the d=1 algorithm: swarm listings (default,
   mothur, uclust), representative seeds, internal structure, per-swarm
   statistics and the amplicon network dump. */

#include "algod1_internal.hpp"
#include <vector>

struct Parameters;  // defined in swarm.hpp
class Data;         // defined in db.hpp


auto write_network_file(uint64_t number_of_networks,
                        struct Parameters const & parameters,
                        Data const & data,
                        std::vector<struct ampinfo_s> const & ampinfo_v,
                        std::vector<unsigned int> & network_v) -> void;

auto output_results(struct Parameters const & parameters,
                    Data const & data,
                    std::vector<struct ampinfo_s> const & ampinfo_v,
                    std::vector<struct swarminfo_s> const & swarminfo_v,
                    Overall_stats const & overall_stats) -> void;

#endif  // SWARM_UTILS_ALGOD1_OUTPUT_H
