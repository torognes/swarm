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

#include "../../utils/fatal.h"
#include "../../utils/system_memory.h"
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // size_t
#include <sys/resource.h>  // Linux: getrusage
// #include <bits/types/struct_rusage.h>  // rusage (since 2017)
#include <sys/sysinfo.h>  // sysinfo
#include <unistd.h>  // sysconf, _SC_PHYS_PAGES, _SC_PAGESIZE


auto system_get_memused() -> uint64_t {
  struct rusage r_usage;  // refactoring: add initializer '{}' (warning with GCC < 5)
  getrusage(RUSAGE_SELF, & r_usage);
  /* Linux: ru_maxrss gives the size in kilobytes  */
  static constexpr unsigned int one_kilobyte {1U << 10U};
  // ru_maxrss is the POSIX-specified 'long' field; glibc nests it in
  // an anonymous union with an internal __ru_maxrss_word alias purely
  // for kernel-ABI flexibility, so this is not a discriminated-union
  // access despite what cppcoreguidelines-pro-type-union-access infers.
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-union-access)
  return static_cast<uint64_t>(r_usage.ru_maxrss * one_kilobyte);
}


auto system_get_memtotal() -> uint64_t {
#if defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  int64_t const phys_pages = sysconf(_SC_PHYS_PAGES);
  int64_t const pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) or (pagesize == -1)) {
    fatal("Cannot determine amount of RAM.");
  }
  return static_cast<uint64_t>(pagesize * phys_pages);

#else

  struct sysinfo info;  // refactoring: add initializer '{}' (warning with GCC < 5)
  if (sysinfo(&info) != 0) {
    fatal("Cannot determine amount of RAM.");
  }
  return info.totalram * info.mem_unit;

#endif
}
