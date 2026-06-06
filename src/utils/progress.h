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

#include <cstdint>  // uint64_t
#include <cstdio>   // std::FILE


struct Parameters;


class Progress
{
public:

  Progress(char const * prompt_, uint64_t size_,
           struct Parameters const & parameters);

  // update() is called from worker threads (algod1.cc network/heavy/light
  // workers) as well as from the main thread; concurrent-safety is
  // guaranteed by each worker holding its own state.mutex when invoking
  // this method.
  auto update(uint64_t current) -> void;

  // Convenience method for single-threaded callers: increments the
  // internal counter by one and reports the new value. Not safe to
  // call from multiple threads concurrently; multi-threaded callers
  // must keep using the explicit-counter update() above.
  auto increment() -> void;

  auto done() const -> void;

private:

  char const * prompt {nullptr};
  uint64_t next {0};
  uint64_t size {0};
  uint64_t chunk {0};
  uint64_t counter {0};
  std::FILE * logfile {nullptr};
  bool silent {false};  // true when output goes to a log file (--log)
};
