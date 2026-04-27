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

#include "swarm.h"
#include "algo.h"
#include "algod1.h"
#include "cli.h"
#include "db.h"
#include "derep.h"
#include "utils/seqinfo.h"
#include "zobrist.h"
#include <cstdint>  // uint64_t
#include <vector>


auto main(int argc, char** argv) -> int
{
  // initialization and checks
  auto const parameters = parse_command_line(argc, argv);

  // parse fasta input
  std::vector<char> data_v;  // refactoring: std::string fails? .data() -> const char *  // alignas(8) does not fix alignment issue
  std::vector<struct seqinfo_s> seqindex_v;
  std::vector<uint64_t> zobrist_tab_base_v;
  std::vector<uint64_t> zobrist_tab_byte_base_v;
  db_read(parameters,
          data_v,
          seqindex_v,
          zobrist_tab_base_v,
          zobrist_tab_byte_base_v);

  // clustering
  switch (parameters.opt_differences)
    {
    case 0:
      dereplicate(parameters);
      break;

    case 1:
      algo_d1_run(parameters);
      break;

    default:
      algo_run(parameters, seqindex_v);
      break;
    }

  // clean up (open output files are closed via RAII when `parameters`
  // goes out of scope)
  db_free();
}
