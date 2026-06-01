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

#ifndef SWARM_UTILS_ALGOD1_STATISTICS_H
#define SWARM_UTILS_ALGOD1_STATISTICS_H


/* Statistics for the d=1 algorithm: classify swarms into light/heavy,
   derive the fastidious Bloom filter geometry, and log the final
   summary. */

#include "algod1_internal.h"
#include <cstdint>  // uint64_t
#include <vector>

struct Parameters;  // defined in swarm.h


auto count_cluster_stats(struct Parameters const & parameters,
                         unsigned int amplicon_count,
                         std::vector<struct swarminfo_s> const & swarminfo_v) -> Cluster_stats;

auto compute_bloom_geometry(struct Parameters const & parameters,
                            uint64_t nucleotides_in_small_clusters) -> Bloom_geometry;

auto log_swarm_summary(struct Parameters const & parameters,
                       Overall_stats const & overall_stats) -> void;

#endif  // SWARM_UTILS_ALGOD1_STATISTICS_H
