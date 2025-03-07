/*
    SWARM

    Copyright (C) 2012-2025 Torbjorn Rognes and Frederic Mahe

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
#include "opt_log.h"
#include "opt_logfile.h"
#include <cstdio>  // fflush, fprintf
#include <cstdint>  // uint64_t


static const char * progress_prompt;
static uint64_t progress_next;
static uint64_t progress_size;
static uint64_t progress_chunk;

auto progress_init(const char * prompt, const uint64_t size) -> void
{
  static constexpr uint64_t progress_granularity {200};
  progress_prompt = prompt;
  progress_size = size;
  progress_chunk = size < progress_granularity ?
    1 : size / progress_granularity;
  progress_next = 1;
  if (not opt_log.empty()) {
    std::fprintf(logfile, "%s", prompt);
  }
  else {
    std::fprintf(logfile, "%s %.0f%%", prompt, 0.0);
  }
}


// refactoring: there are three calls to this function that are beyond
// the pthread wall. There is no easy way (for now) to pass additional
// arguments beyond that wall. This is a major roadblock and prevents
// us to eliminate global variables (opt_log, logfile, as well as
// 'progress_*' global). Could be solved by using std::thread?
auto progress_update(const uint64_t progress) -> void
{
  if (opt_log.empty() and (progress >= progress_next))
    {
      std::fprintf(logfile, "  \r%s %.0f%%", progress_prompt,
                   100.0 * static_cast<double>(progress)
                   / static_cast<double>(progress_size));
      progress_next = progress + progress_chunk;
      std::fflush(logfile);
    }
}


auto progress_done(struct Parameters const & parameters) -> void
{
  if (not parameters.opt_log.empty()) {
    std::fprintf(parameters.logfile, " %.0f%%\n", 100.0);
  }
  else {
    std::fprintf(parameters.logfile, "  \r%s %.0f%%\n", progress_prompt, 100.0);
  }
  std::fflush(parameters.logfile);
}
