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

#include "../../utils/fatal.hpp"
#include "../../utils/system_memory.hpp"
#include <cstdint>  // uint64_t
#include <windows.h>
#include <psapi.h>


auto system_get_memused() -> uint64_t {
  PROCESS_MEMORY_COUNTERS pmc {};
  // a zero return means failure; bail out as the POSIX backends do
  // rather than returning the uninitialised struct's garbage.
  if (GetProcessMemoryInfo(GetCurrentProcess(),
                           &pmc,
                           sizeof(PROCESS_MEMORY_COUNTERS)) == 0) {
    fatal("Cannot determine amount of RAM.");
  }
  return pmc.PeakWorkingSetSize;
}


auto system_get_memtotal() -> uint64_t {
  MEMORYSTATUSEX memory_status {};
  memory_status.dwLength = sizeof(MEMORYSTATUSEX);
  if (GlobalMemoryStatusEx(&memory_status) == 0) {
    fatal("Cannot determine amount of RAM.");
  }
  return memory_status.ullTotalPhys;
}
