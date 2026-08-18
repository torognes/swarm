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

#ifndef SWARM_UTILS_SEQ_INDEX_H
#define SWARM_UTILS_SEQ_INDEX_H


#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t


// (offset, length) pair into a flat buffer; used to locate headers
// and packed sequences inside the two database vectors. The unit is
// the buffer's element: for a header, offset counts bytes into the
// header vector and length is the byte count; for a sequence, offset
// counts 64-bit words into the packed-sequence vector and length is
// the nucleotide count.
struct Index {
  uint64_t offset {0};
  std::size_t length {0};
};


// One fasta record as captured during parsing: source line number
// plus header and sequence locations in the database vectors.
struct Entry {
  // 64-bit: the line count is not bounded by the sequence count -- one
  // long wrapped record spans many lines -- and everything downstream
  // (find_abundance, Seq_stats::missingabundance_lineno) is already
  // uint64_t. Widening costs nothing here: the two Index members align
  // the struct to 8 bytes, so a 32-bit lineno paid the same 8 bytes in
  // padding.
  uint64_t lineno {1};
  struct Index header;
  struct Index sequence;
};

#endif  // SWARM_UTILS_SEQ_INDEX_H
