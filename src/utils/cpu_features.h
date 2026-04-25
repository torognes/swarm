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

#ifndef SWARM_UTILS_CPU_FEATURES_H
#define SWARM_UTILS_CPU_FEATURES_H

// Subset of CPU feature flags consulted from hot paths (search8/search16
// dispatch on ssse3/sse41, qgram comparison dispatches on popcnt). Bundled
// here to avoid pulling the full Parameters struct into low-level files.
//
// Aggregate type (no default member initializers, no constructor) so it can
// be brace-initialized in C++11.
struct Cpu_features
{
  bool ssse3;
  bool sse41;
  bool popcnt;
};

#endif
