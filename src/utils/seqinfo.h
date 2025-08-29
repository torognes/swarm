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

#include <cstdint>  // uint64_t

// refactoring: add View header_view and View seq_view to
// progressively get rid of raw pointers

struct seqinfo_s
{
  char const * header;
  char const * seq;
  uint64_t abundance;
  uint64_t hdrhash;
  uint64_t seqhash;
  int headerlen;  // refactoring: should be unsigned int or size_t
  unsigned int seqlen;
  unsigned int clusterid;
  int abundance_start;
  int abundance_end;
  int dummy; /* alignment padding only */
};
