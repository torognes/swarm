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


auto progress_init(struct Progress_status & progress,
                   char const * prompt, uint64_t const size,
                   struct Parameters const & parameters) -> void
{
  static constexpr uint64_t progress_granularity {200};
  progress.prompt = prompt;
  progress.size = size;
  progress.chunk = size < progress_granularity ?
    1 : size / progress_granularity;
  progress.next = 1;
  progress.logfile = parameters.logfile;
  progress.silent = not parameters.opt_log.empty();
  if (progress.silent) {
    std::fprintf(progress.logfile, "%s", prompt);
  }
  else {
    std::fprintf(progress.logfile, "%s %.0f%%", prompt, 0.0);
  }
}


// Called from within worker threads (algod1.cc network/heavy/light workers)
// as well as from the main thread; concurrent-safety is guaranteed by each
// worker holding its own state.mutex when invoking this function.
auto progress_update(struct Progress_status & progress, uint64_t const current) -> void
{
  if (progress.silent) { return; }  // no progress output if log is a file
  if (current < progress.next) { return; }  // milestone not yet reached
  std::fprintf(progress.logfile, "  \r%s %.0f%%", progress.prompt,
               100.0 * static_cast<double>(current)
               / static_cast<double>(progress.size));
  progress.next = current + progress.chunk;
  std::fflush(progress.logfile);
}


auto progress_done(struct Progress_status const & progress) -> void
{
  if (progress.silent) {
    std::fprintf(progress.logfile, " %.0f%%\n", 100.0);
  }
  else {
    std::fprintf(progress.logfile, "  \r%s %.0f%%\n", progress.prompt, 100.0);
  }
  std::fflush(progress.logfile);
}
