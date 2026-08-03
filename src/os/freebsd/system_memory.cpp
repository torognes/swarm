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
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <sys/resource.h>  // getrusage, RUSAGE_SELF, struct rusage
#include <sys/sysctl.h>  // sysctlbyname


auto system_get_memused() -> uint64_t {
  struct rusage r_usage;  // refactoring: add initializer '{}' (warning with GCC < 5)
  getrusage(RUSAGE_SELF, & r_usage);
  /* FreeBSD: ru_maxrss is in kilobytes (same as Linux) */
  static constexpr unsigned int one_kilobyte {1U << 10U};
  // ru_maxrss is the POSIX-specified 'long' field; libc may nest it
  // in an anonymous union for kernel-ABI flexibility, so this is not
  // a discriminated-union access despite what
  // cppcoreguidelines-pro-type-union-access infers.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-union-access)
  return static_cast<uint64_t>(r_usage.ru_maxrss * one_kilobyte);
}


auto system_get_memtotal() -> uint64_t {
  /* sysctlbyname("hw.physmem") writes a uint64_t directly, avoiding
     the 32-bit overflow you would get from the older
     sysctl({CTL_HW, HW_PHYSMEM}) interface on hosts with >= 4 GB. */
  uint64_t ram {0};
  std::size_t length = sizeof(ram);
  if (sysctlbyname("hw.physmem", &ram, &length, nullptr, 0) != 0) {
    fatal("Cannot determine amount of RAM.");
  }
  return ram;
}
