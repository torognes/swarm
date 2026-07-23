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

#include "view.hpp"
#include <cstdint>  // uint64_t

// refactoring: header (char const *) + headerlen (int) merged into
// header_view (View<char>); seq + seqlen could be similarly
// merged in a follow-up, but the byte/nucleotide-count mismatch
// (4 nt packed per byte) makes a clean swap less obvious.

struct seqinfo_s
{
  View<char> header_view;
  char const * seq {nullptr};
  uint64_t abundance {0};
  uint64_t seqhash {0};
  unsigned int seqlen {0};
  int abundance_start {0};
  int abundance_end {0};
};
