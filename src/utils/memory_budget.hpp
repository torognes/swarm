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

#ifndef SWARM_UTILS_MEMORY_BUDGET_H
#define SWARM_UTILS_MEMORY_BUDGET_H

#include "fatal.hpp"
#include "system_memory.hpp"
#include <cstdint>  // uint64_t


// Guard an allocation that scales super-linearly with the input (the
// O(L^2) alignment matrices used when d > 1). Swarm builds with
// -fno-exceptions, so a std::vector allocation failure calls
// std::terminate() rather than throwing; this converts the pathological
// case (a very long input sequence) into an actionable error *before* the
// request reaches operator new.
//
// The check is expressed as a division so 'bytes_per_unit * units' cannot
// overflow uint64_t for a maximum-length sequence.
inline auto require_ram(uint64_t const bytes_per_unit,
                        uint64_t const units,
                        char const * const context) -> void {
  if (units == 0) { return; }
  if (bytes_per_unit > system_get_memtotal() / units) {
    static constexpr uint64_t bytes_per_mib {uint64_t{1} << 20U};
    fatal("Not enough memory for ", context, ": about ",
          (bytes_per_unit / bytes_per_mib) * units,
          " MB is required, which exceeds the total amount of RAM. Reduce "
          "the number of threads (-t) and/or the length of the input "
          "sequences.");
  }
}

#endif  // SWARM_UTILS_MEMORY_BUDGET_H
