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

#include "../swarm.h"
#include "progress.h"
#include <cstdio>  // fflush, fprintf
#include <cstdint>  // uint64_t


namespace {
  constexpr uint64_t progress_granularity {200};
}


Progress::Progress(char const * prompt_, uint64_t const size_,
                   struct Parameters const & parameters)
  : prompt(prompt_),
    next(1),
    size(size_),
    chunk(size_ < progress_granularity ? 1 : size_ / progress_granularity),
    logfile(parameters.logfile),
    silent(not parameters.opt_log.empty()) {
  if (silent) {
    std::fprintf(logfile, "%s", prompt);
  }
  else {
    std::fprintf(logfile, "%s %.0f%%", prompt, 0.0);
  }
}


auto Progress::update(uint64_t const current) -> void {
  if (silent) { return; }  // no progress output if log is a file
  if (current < next) { return; }  // milestone not yet reached
  std::fprintf(logfile, "  \r%s %.0f%%", prompt,
               100.0 * static_cast<double>(current)
               / static_cast<double>(size));
  next = current + chunk;
  std::fflush(logfile);
}


auto Progress::increment() -> void {
  ++counter;
  update(counter);
}


auto Progress::done() const -> void {
  if (silent) {
    std::fprintf(logfile, " %.0f%%\n", 100.0);
  }
  else {
    std::fprintf(logfile, "  \r%s %.0f%%\n", prompt, 100.0);
  }
  std::fflush(logfile);
}
