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

#include "../swarm.hpp"
#include "progress.hpp"
#include <cstdio>  // fflush, fprintf, fputs
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
  // Both branches open with the prompt; only the percentage differs, and
  // "%.0f" of the literal 0.0 is always the one character '0', so there is
  // nothing here for a format string to decide.
  static_cast<void>(std::fputs(prompt, logfile));
  if (not silent) {
    static_cast<void>(std::fputs(" 0%", logfile));
  }
}


auto Progress::update(uint64_t const current) -> void {
  if (silent) { return; }  // no progress output if log is a file
  if (current < next) { return; }  // milestone not yet reached
  // The one fprintf swarm still needs: "%.0f" rounds half to even, which
  // integer truncation does not reproduce (1.5 % -> "2" against "1"), and
  // the milestones land on half-integer percentages because chunk is
  // size/200. Matching it exactly would take more code than the double.
  //
  // LTO inlines update() into ~29 callers, so this fprintf is duplicated
  // that many times. Measured: hiding it behind a noinline helper cuts the
  // copies from 33 to 5 and .text by 23 bytes. Not worth a function and a
  // compiler-specific attribute; recorded so it is not re-attempted.
  static_cast<void>(std::fprintf(logfile, "  \r%s %.0f%%", prompt,
                                 100.0 * static_cast<double>(current)
                                 / static_cast<double>(size)));
  next = current + chunk;
  std::fflush(logfile);
}


auto Progress::increment() -> void {
  ++counter;
  update(counter);
}


auto Progress::done() const -> void {
  // Same shape as the constructor: the percentage is the literal 100.0, so
  // "%.0f" always yields "100". The non-silent branch additionally rewinds
  // over the partial line update() left behind.
  if (not silent) {
    static_cast<void>(std::fputs("  \r", logfile));
    static_cast<void>(std::fputs(prompt, logfile));
  }
  static_cast<void>(std::fputs(" 100%\n", logfile));
  std::fflush(logfile);
}
