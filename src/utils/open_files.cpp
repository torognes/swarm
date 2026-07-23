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


auto open_files(struct Parameters & parameters) -> void {
  // special case (always '-')??
  parameters.outfile = fopen_output(parameters.opt_output_file.c_str());
  if (not parameters.outfile) {
    fatal("Unable to open output file for writing.");
  }

  /* open files */

  if (not parameters.opt_log.empty()) {
      parameters.logfile_handle = fopen_output(parameters.opt_log.c_str());
      if (not parameters.logfile_handle) {
        fatal("Unable to open log file for writing.");
      }
      parameters.logfile = parameters.logfile_handle.get();
  }

  if (not parameters.opt_seeds.empty()) {
      parameters.seeds_file = fopen_output(parameters.opt_seeds.c_str());
      if (not parameters.seeds_file) {
        fatal("Unable to open seeds file for writing.");
      }
  }

  if (not parameters.opt_statistics_file.empty()) {
      parameters.statsfile = fopen_output(parameters.opt_statistics_file.c_str());
      if (not parameters.statsfile) {
        fatal("Unable to open statistics file for writing.");
      }
  }

  if (not parameters.opt_uclust_file.empty()) {
      parameters.uclustfile = fopen_output(parameters.opt_uclust_file.c_str());
      if (not parameters.uclustfile) {
        fatal("Unable to open uclust file for writing.");
      }
  }

  if (not parameters.opt_internal_structure.empty()) {
      parameters.internal_structure_file = fopen_output(parameters.opt_internal_structure.c_str());
      if (not parameters.internal_structure_file) {
        fatal("Unable to open internal structure file for writing.");
      }
  }

  if (not parameters.opt_network_file.empty()) {
      parameters.network_file = fopen_output(parameters.opt_network_file.c_str());
      if (not parameters.network_file) {
        fatal("Unable to open network file for writing.");
      }
  }
}
