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

#include "open_files.hpp"
#include "../swarm.hpp"
#include "fatal.hpp"
#include "input_output.hpp"

#include <array>  // std::array
#include <string>  // std::string


namespace {

  // One optional output file: the option that names it, the handle that
  // owns it, and what to say when it cannot be opened. The option string
  // and the handle are references into Parameters, so the loop below
  // writes the opened handle back where the rest of swarm reads it. A
  // reference member does not inherit the constness of the element it sits
  // in, which is why the table itself can be const.
  struct Optional_output_file {
    std::string const & option;  // empty unless the option was given
    FileHandle & handle;
    char const * message;
  };

}  // end of anonymous namespace


auto open_files(struct Parameters & parameters) -> void {
  // Not part of the table below, and not conditional: opt_output_file
  // defaults to '-', so the main output is always opened, and an empty name
  // must reach fopen_output rather than being skipped as "not requested"
  // (cli.cpp rejects it before that, so what remains here is a name fopen
  // refused).
  parameters.outfile = fopen_output(parameters.opt_output_file);
  if (not parameters.outfile) {
    fatal("Unable to open output file for writing.");
  }

  /* open files */

  std::array<Optional_output_file, 6> const optional_files {{
      {parameters.opt_log, parameters.logfile_handle,
       "Unable to open log file for writing."},
      {parameters.opt_seeds, parameters.seeds_file,
       "Unable to open seeds file for writing."},
      {parameters.opt_statistics_file, parameters.statsfile,
       "Unable to open statistics file for writing."},
      {parameters.opt_uclust_file, parameters.uclustfile,
       "Unable to open uclust file for writing."},
      {parameters.opt_internal_structure, parameters.internal_structure_file,
       "Unable to open internal structure file for writing."},
      {parameters.opt_network_file, parameters.network_file,
       "Unable to open network file for writing."}
    }};

  for (auto const & entry : optional_files) {
    if (entry.option.empty()) { continue; }
    entry.handle = fopen_output(entry.option);
    if (not entry.handle) {
      fatal(entry.message);
    }
  }

  // logfile is a raw alias, not an owner: it defaults to stderr and only
  // points into logfile_handle when -l asked for a file.
  if (parameters.logfile_handle) {
    parameters.logfile = parameters.logfile_handle.get();
  }
}
