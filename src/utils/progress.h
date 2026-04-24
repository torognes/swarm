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


struct Progress_status
{
  char const * prompt {nullptr};
  uint64_t next {0};
  uint64_t size {0};
  uint64_t chunk {0};
  std::FILE * logfile {nullptr};
  bool silent {false};  // true when output goes to a log file (--log)
};


auto progress_init(struct Progress_status & progress,
                   char const * prompt, uint64_t size,
                   struct Parameters const & parameters) -> void;
auto progress_update(struct Progress_status & progress, uint64_t current) -> void;
auto progress_done(struct Progress_status const & progress) -> void;
