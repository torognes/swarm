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

#ifndef SWARM_UTILS_SYSTEM_MEMORY_H
#define SWARM_UTILS_SYSTEM_MEMORY_H


#include <cstdint>


// operating system specific functions (Windows, macOS and Linux)
auto system_get_memused() -> uint64_t;

// The machine's memory. Prefer system_get_memlimit() below for anything
// that budgets: this figure is the host's, and a container's cgroup limit is
// invisible in it, so inside one it can be orders of magnitude too large.
auto system_get_memtotal() -> uint64_t;

// The memory this process may actually use: system_get_memtotal(), reduced
// to the smallest limit that applies to it. On Linux that is the cgroup
// memory limit, which is how Slurm, Docker, podman and Kubernetes cap a job
// -- they share the host's kernel, so there is no smaller "total" for
// system_get_memtotal() to report. Platforms with no such mechanism visible
// return system_get_memtotal() unchanged.
//
// Note that a virtual machine needs nothing here: it has its own kernel, and
// system_get_memtotal() already reports what the hypervisor gave the guest.
//
// Never zero, and never larger than system_get_memtotal().
auto system_get_memlimit() -> uint64_t;

#endif  // SWARM_UTILS_SYSTEM_MEMORY_H
