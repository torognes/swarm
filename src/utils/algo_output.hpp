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

#ifndef SWARM_UTILS_ALGO_OUTPUT_H
#define SWARM_UTILS_ALGO_OUTPUT_H


/* Result writers for the d > 1 algorithm: swarm listings (default,
   mothur), the per-amplicon internal structure, representative seeds in
   fasta format, and the per-swarm uclust and statistics records. */

#include "algo_internal.hpp"  // ampliconinfo_s, Cluster_state
#include "view.hpp"         // View
#include <cstdint>          // uint64_t

struct Parameters;  // defined in swarm.hpp
class Data;         // defined in db.hpp
class NwAligner;    // defined in utils/nw_aligner.hpp


// amps is the whole amplicon pool, in cluster order. It carries its own
// length, so there is no separate amplicon count that a later edit could
// leave disagreeing with it.
auto write_swarms_default_format(struct Parameters const & parameters,
                                 Data const & data,
                                 View<struct ampliconinfo_s> amps) -> void;

auto write_swarms_mothur_format(unsigned int swarmid,
                                struct Parameters const & parameters,
                                Data const & data,
                                View<struct ampliconinfo_s> amps) -> void;

auto write_internal_structure_line(uint64_t parent_id,
                                   uint64_t child_id,
                                   uint64_t diff,
                                   unsigned int swarmid,
                                   unsigned int generation,
                                   struct Parameters const & parameters,
                                   Data const & data) -> void;

auto write_representative_sequences(struct Parameters const & parameters,
                                    Data const & data,
                                    View<struct ampliconinfo_s> amps) -> void;

// hits is the cluster's members, seed first: the filled part of the
// caller's hit buffer, not the pool-sized buffer behind it.
auto write_cluster_outputs(unsigned int swarmid,
                           uint64_t seedampliconid,
                           Cluster_state const & state,
                           View<uint64_t> hits,
                           NwAligner * aligner,
                           struct Parameters const & parameters,
                           Data const & data) -> void;

#endif  // SWARM_UTILS_ALGO_OUTPUT_H
