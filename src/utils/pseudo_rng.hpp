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

#ifndef SWARM_UTILS_PSEUDO_RNG_H
#define SWARM_UTILS_PSEUDO_RNG_H


#include <cstdint>
#include <random>

// pseudo random number generator:
// Mersenne Twister uint64 uniform distribution
// (initialized only once for reproducibility,
//  then each call produces a distinct uint64 value)
//
// The engine lives in a function-local static, which
//   - sidesteps cppcoreguidelines-avoid-non-const-global-variables
//     (the rule only flags namespace-scope mutable variables),
//   - preserves the previous "one independent seed=1 stream per
//     translation unit" semantic: 'static' on the function gives it
//     internal linkage, so every TU that includes this header keeps
//     its own engine state.
constexpr static unsigned int seed {1};
static auto rand_64() -> uint64_t {
  static std::mt19937_64 engine(seed);
  return engine();
}

#endif  // SWARM_UTILS_PSEUDO_RNG_H
