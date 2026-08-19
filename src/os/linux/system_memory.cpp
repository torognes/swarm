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
#include <algorithm>  // std::min
#include <array>
#include <cerrno>  // errno, ERANGE
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // std::FILE, std::fopen, std::fgets, std::fclose
#include <cstdlib>  // std::strtoull
#include <memory>  // std::unique_ptr
#include <string>
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


// The cgroup memory limit, which is what actually caps a containerised run
// (Slurm with cgroups, Docker, podman, Kubernetes, Apptainer). A container
// shares the host's kernel, so sysconf(_SC_PHYS_PAGES) above reports the
// host's memory and knows nothing about the limit; a virtual machine, by
// contrast, has its own kernel and needs none of this.
namespace {

  constexpr uint64_t no_limit {0};

  // cgroup v1 has no word for "unlimited": it writes a sentinel close to
  // INT64_MAX rounded down to a page multiple (9223372036854771712). Read
  // literally that is an 8 EiB budget, which would be worse than no cap at
  // all, so anything this large means "no limit".
  constexpr uint64_t implausible_limit {uint64_t{1} << 62U};

  // RAII for the two /proc and /sys reads below. Not input_output.hpp's
  // FileHandle: that one's deleter calls fatal() when the stream carries an
  // error, which is the right policy for swarm's own output files and the
  // wrong one for an optional probe of a kernel interface.
  struct Close_file {
    auto operator()(std::FILE * const handle) const -> void {
      static_cast<void>(std::fclose(handle));
    }
  };
  using Probe_file = std::unique_ptr<std::FILE, Close_file>;

  // One line, however long: cgroup paths under Kubernetes run to a couple of
  // hundred characters, and a truncated path would silently name a file that
  // does not exist. Returns false at end of input.
  auto read_line(std::FILE * const input, std::string & line) -> bool {
    line.clear();
    std::array<char, 256> buffer {};
    while (std::fgets(buffer.data(), static_cast<int>(buffer.size()), input) != nullptr) {
      line.append(buffer.data());
      if (not line.empty() and (line.back() == '\n')) {
        line.pop_back();
        return true;
      }
    }
    return not line.empty();
  }

  // The limit in one cgroup file, or no_limit when the file is absent, holds
  // cgroup v2's literal "max", or cannot be parsed.
  auto read_limit_file(std::string const & path) -> uint64_t {
    Probe_file const input {std::fopen(path.c_str(), "r")};
    if (not input) { return no_limit; }

    std::string line;
    if (not read_line(input.get(), line)) { return no_limit; }

    static constexpr int base_value {10};
    char * end_ptr {nullptr};
    errno = 0;
    auto const value = std::strtoull(line.c_str(), &end_ptr, base_value);
    if (errno == ERANGE) { return no_limit; }
    if (end_ptr == line.c_str()) { return no_limit; }  // "max", or not a number
    if (*end_ptr != '\0') { return no_limit; }          // trailing junk
    if (value >= implausible_limit) { return no_limit; }
    return value;
  }

  // Where this process sits in the hierarchy, from /proc/self/cgroup. The
  // unified hierarchy (cgroup v2) writes a single "0::/path" line; cgroup v1
  // writes one line per controller as "id:controllers:/path", and the one
  // naming the memory controller is the relevant one.
  auto own_cgroup_path(char const * const proc_file, bool & unified) -> std::string {
    unified = false;
    Probe_file const input {std::fopen(proc_file, "r")};
    if (not input) { return std::string(); }

    std::string line;
    std::string memory_path;
    while (read_line(input.get(), line)) {
      if (line.compare(0, 3, "0::") == 0) {
        unified = true;
        return line.substr(3);
      }
      auto const first = line.find(':');
      if (first == std::string::npos) { continue; }
      auto const second = line.find(':', first + 1);
      if (second == std::string::npos) { continue; }
      if (line.find("memory", first + 1) < second) {
        memory_path = line.substr(second + 1);
      }
    }
    return memory_path;
  }

  // The smallest limit set anywhere on this process's cgroup path.
  //
  // The walk *upward* is the point, not a refinement: Kubernetes sets the
  // limit on the pod's cgroup, an ancestor of the container's, whose own
  // memory.max reads "max". Reading the leaf alone finds nothing and
  // concludes there is no limit, which is the bug this function exists to
  // fix. The tightest limit on the path is the one that binds.
  //
  // The two paths are parameters so that the parsing and the walk can be
  // exercised against a synthetic hierarchy: neither swarm-tests nor an
  // unprivileged CI job can create a real cgroup, and the case a real host
  // usually presents is the one that returns no_limit.
  auto cgroup_memory_limit(std::string const & cgroup_root,
                           char const * const proc_file) -> uint64_t {
    bool unified {false};
    auto path = own_cgroup_path(proc_file, unified);
    auto const base = unified ? cgroup_root : cgroup_root + "/memory";
    std::string const leaf = unified ? "/memory.max" : "/memory.limit_in_bytes";

    auto limit = no_limit;
    while (true) {
      auto const value = read_limit_file(base + path + leaf);
      if ((value != no_limit) and ((limit == no_limit) or (value < limit))) {
        limit = value;
      }
      auto const slash = path.rfind('/');
      if (slash == std::string::npos) { break; }
      path.erase(slash);
    }
    return limit;
  }

}  // namespace


auto system_get_memlimit() -> uint64_t {
  auto const physical = system_get_memtotal();
  auto const limit = cgroup_memory_limit("/sys/fs/cgroup", "/proc/self/cgroup");
  if (limit == no_limit) { return physical; }
  return std::min(limit, physical);
}
