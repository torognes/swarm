/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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

#include "utils/fatal.h"
#include <cstdint>  // int64_t, uint64_t
#include <cstdio> // size_t

#ifdef __APPLE__
#include <sys/resource.h>
#include <sys/sysctl.h>
#elif defined _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>  // Linux: getrusage
#include <bits/types/struct_rusage.h>  // rusage (since 2017)
#include <sys/sysinfo.h>  // sysinfo
#endif


auto arch_get_memused() -> uint64_t
{
#ifdef _WIN32

  PROCESS_MEMORY_COUNTERS pmc;
  GetProcessMemoryInfo(GetCurrentProcess(),
                       &pmc,
                       sizeof(PROCESS_MEMORY_COUNTERS));
  return pmc.PeakWorkingSetSize;

#else

  struct rusage r_usage {};
  getrusage(RUSAGE_SELF, & r_usage);

# ifdef __APPLE__
  /* Mac: ru_maxrss gives the size in bytes */
  return static_cast<uint64_t>(r_usage.ru_maxrss);
# else
  /* Linux: ru_maxrss gives the size in kilobytes  */
  static constexpr unsigned int one_kilobyte {1U << 10U};
  return static_cast<uint64_t>(r_usage.ru_maxrss * one_kilobyte);
# endif

#endif
}


auto arch_get_memtotal() -> uint64_t
{
#ifdef _WIN32

  MEMORYSTATUSEX memory_status;
  memory_status.dwLength = sizeof(MEMORYSTATUSEX);
  GlobalMemoryStatusEx(&memory_status);
  return memory_status.ullTotalPhys;

#elif defined(__APPLE__)

  int mib [] = { CTL_HW, HW_MEMSIZE };
  int64_t ram = 0;
  std::size_t length = sizeof(ram);
  if(sysctl(mib, 2, &ram, &length, nullptr, 0) != 0) {
    fatal(error_prefix, "Cannot determine amount of RAM.");
  }
  return static_cast<uint64_t>(ram);

#elif defined(_SC_PHYS_PAGES) && defined(_SC_PAGESIZE)

  const int64_t phys_pages = sysconf(_SC_PHYS_PAGES);
  const int64_t pagesize = sysconf(_SC_PAGESIZE);
  if ((phys_pages == -1) || (pagesize == -1)) {
    fatal(error_prefix, "Cannot determine amount of RAM.");
  }
  return static_cast<uint64_t>(pagesize * phys_pages);

#else

  struct sysinfo info {};
  if (sysinfo(&info) != 0) {
    fatal(error_prefix, "Cannot determine amount of RAM.");
  }
  return info.totalram * info.mem_unit;

#endif
}
