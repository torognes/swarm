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
#include "fatal.h"
#include "input_output.h"
#include "opt_logfile.h"
#include <cstdio>  // FILE, fclose, stderr  // refactoring: replace with <fstream>
#include <vector>


std::FILE * logfile {stderr};  // cstdio stderr macro is expanded to type std::FILE*


auto open_files(struct Parameters & parameters) -> void
{
  // special case (always '-')??
  parameters.outfile = fopen_output(parameters.opt_output_file.c_str());
  if (parameters.outfile == nullptr) {
    fatal(error_prefix, "Unable to open output file for writing.");
  }

  /* open files */

  if (not parameters.opt_log.empty())
    {
      parameters.logfile = fopen_output(parameters.opt_log.c_str());
      logfile = parameters.logfile;
      if (parameters.logfile == nullptr) {
        fatal(error_prefix, "Unable to open log file for writing.");
      }
    }

  if (not parameters.opt_seeds.empty())
    {
      parameters.seeds_file = fopen_output(parameters.opt_seeds.c_str());
      if (parameters.seeds_file == nullptr) {
        fatal(error_prefix, "Unable to open seeds file for writing.");
      }
    }

  if (not parameters.opt_statistics_file.empty())
    {
      parameters.statsfile = fopen_output(parameters.opt_statistics_file.c_str());
      if (parameters.statsfile == nullptr) {
        fatal(error_prefix, "Unable to open statistics file for writing.");
      }
    }

  if (not parameters.opt_uclust_file.empty())
    {
      parameters.uclustfile = fopen_output(parameters.opt_uclust_file.c_str());
      if (parameters.uclustfile == nullptr) {
        fatal(error_prefix, "Unable to open uclust file for writing.");
      }
    }

  if (not parameters.opt_internal_structure.empty())
    {
      parameters.internal_structure_file = fopen_output(parameters.opt_internal_structure.c_str());
      if (parameters.internal_structure_file == nullptr) {
        fatal(error_prefix, "Unable to open internal structure file for writing.");
      }
    }

  if (not parameters.opt_network_file.empty())
    {
      parameters.network_file = fopen_output(parameters.opt_network_file.c_str());
      if (parameters.network_file == nullptr) {
        fatal(error_prefix, "Unable to open network file for writing.");
      }
    }
}


auto close_files(struct Parameters & parameters) -> void {
  const std::vector<std::FILE *> file_handles
    {parameters.network_file, parameters.internal_structure_file,
     parameters.uclustfile, parameters.statsfile, parameters.seeds_file, parameters.outfile,
     parameters.logfile};
  for (auto * const file_handle : file_handles) {
    if (file_handle != nullptr) {
      std::fclose(file_handle);
    }
  }
}
