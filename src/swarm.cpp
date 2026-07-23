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

#include "swarm.hpp"  // struct Parameters
#include "algo.hpp"
#include "algod1.hpp"
#include "cli.hpp"  // parse_command_line()
#include "db.hpp"  // class Data
#include "dereplicate.hpp"


auto main(int const argc, char * const * argv) -> int {
  // initialization and checks
  auto const parameters = parse_command_line(argc, argv);

  // parse fasta input
  Data const data {parameters};

  // clustering
  switch (parameters.opt_differences) {
  case 0:
    dereplicate(parameters, data);
    break;

  case 1:
    algo_d1_run(parameters, data);
    break;

  default:
    algo_run(parameters, data);
    break;
  }
}
