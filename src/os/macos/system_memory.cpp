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
#include <array>  // std::array
#include <cstddef>  // std::size_t
#include <cstdint>  // int64_t, uint64_t
#include <sys/resource.h>
#include <sys/sysctl.h>


auto system_get_memused() -> uint64_t {
  struct rusage r_usage;  // refactoring: add initializer '{}' (warning with GCC < 5)
  getrusage(RUSAGE_SELF, & r_usage);
  /* Mac: ru_maxrss gives the size in bytes */
  return static_cast<uint64_t>(r_usage.ru_maxrss);
}


auto system_get_memtotal() -> uint64_t {
  std::array<int, 2> mib {{ CTL_HW, HW_MEMSIZE }};
  int64_t ram = 0;
  std::size_t length = sizeof(ram);
  if (sysctl(mib.data(), static_cast<unsigned int>(mib.size()),
             &ram, &length, nullptr, 0) != 0) {
    fatal("Cannot determine amount of RAM.");
  }
  return static_cast<uint64_t>(ram);
}


// macOS has no cgroup-like per-process memory limit to consult, so the
// machine's memory is also the limit. (Containers on macOS run inside a
// Linux virtual machine, whose guest kernel reports its own assignment
// through the Linux implementation of system_get_memtotal.)
auto system_get_memlimit() -> uint64_t {
  return system_get_memtotal();
}
