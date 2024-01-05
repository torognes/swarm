/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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

#include "util.h"
#include "utils/nt_codec.h"
#include <algorithm>  // std::min()
#include <array>
#include <cassert>  // assert()
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // snprintf, size_t
#include <cstring>  // memcpy, memmove, memset
#include <vector>


constexpr auto n_cells = 32ULL;  // number of chars in sym_nt
constexpr unsigned char maskup      = 1;
constexpr unsigned char maskleft    = 2;
constexpr unsigned char maskextup   = 4;
constexpr unsigned char maskextleft = 8;


auto pushop(const char newop, char * & cigarendp, char & operation, int & count) -> void
{
  static constexpr auto buffer_length = 25U;

  if (newop == operation) {
    ++count;
  }
  else
  {
    --cigarendp;
    *cigarendp = operation;  // write char operation at position end - 1
    if (count > 1)
    {
      std::array<char, buffer_length> buf {{}};
      const auto len = std::snprintf(buf.data(), buffer_length, "%d", count);
      assert(len >= 0);
      cigarendp -= len;
      std::memcpy(cigarendp, buf.data(), static_cast<std::size_t>(len));
    }
    operation = newop;
    count = 1;
  }
}


auto finishop(char * & cigarendp, char & operation, int & count) -> void
{
  static constexpr auto buffer_length = 25U;

  if ((operation != '\0') and (count != 0))
  {
    --cigarendp;
    *cigarendp = operation;  // write char operation at position end - 1
    if (count > 1)
    {
      std::array<char, buffer_length> buf {{}};
      const auto len = std::snprintf(buf.data(), buffer_length, "%d", count);
      assert(len >= 0);
      cigarendp -= len;
      std::memcpy(cigarendp, buf.data(), static_cast<std::size_t>(len));
    }
    operation = '\0';
    count = 0;
  }
}


auto align(char * dseq,
           const uint64_t dlen,
           char * qseq,
           const uint64_t qlen,
           const std::array<int64_t, n_cells * n_cells> & score_matrix,
           const uint64_t gapopen,
           const uint64_t gapextend,
           std::vector<unsigned char> & directions,
           std::vector<uint64_t> & hearray) -> void
{
  // alignment priority when backtracking (from lower right corner):
  // 1. left/insert/e (gap in query sequence (qseq))
  // 2. diagonal/align/h (match/mismatch)
  // 3. top (up)/delete/f (gap in database sequence (dseq))
  static constexpr auto multiplier = 5U;

  assert(directions.size() >= qlen * dlen);
  assert(hearray.size() >= 2 * qlen);

  // hearray: array of alignments and insertions (initialized here,
  // then modified through a pointer; can't be const)
  for(auto column = 0UL; column < qlen; column++)
    {
      hearray[2 * column]     = 1 * gapopen + (column + 1) * gapextend; // H (N)
      hearray[2 * column + 1] = 2 * gapopen + (column + 2) * gapextend; // E
    }

  uint64_t previous_diagonal {0};
  uint64_t left {0};

  for(auto row = 0UL; row < dlen; row++)
    {
      auto he_index = 0UL;
      auto top = 2 * gapopen + (row + 2) * gapextend;
      uint64_t diagonal = (row == 0) ? 0 : (gapopen + row * gapextend);

      for(auto column = 0UL; column < qlen; column++)
        {
          const auto index = qlen * row + column;

          previous_diagonal = hearray[he_index];
          left = hearray[he_index + 1];
          diagonal += static_cast<uint64_t>(score_matrix
                                            [((nt_extract(dseq, row) + 1U) << multiplier)
                                             + (nt_extract(qseq, column) + 1)]);

          directions[index] |= (top < diagonal ? maskup : 0U);
          diagonal = std::min(diagonal, top);
          diagonal = std::min(diagonal, left);
          directions[index] |= (left == diagonal ? maskleft : 0U);

          hearray[he_index] = diagonal;

          diagonal += gapopen + gapextend;
          left += gapextend;
          top += gapextend;

          directions[index] |= (top < diagonal ? maskextup : 0U);
          directions[index] |= (left < diagonal ? maskextleft : 0U);
          top = std::min(diagonal, top);
          left = std::min(diagonal, left);

          hearray[he_index + 1] = left;
          diagonal = previous_diagonal;
          he_index += 2;
        }
    }
}


auto backtrack(char * dseq,
               const uint64_t dlen,
               char * qseq,
               const uint64_t qlen,
               uint64_t & nwdiff,
               uint64_t & nwalignmentlength,
               char * & nwalignment,
               std::vector<unsigned char> const & directions) -> void
{
  /* backtrack: count differences and save alignment in cigar string */

  uint64_t alength {0};
  uint64_t matches {0};

  // refactoring: complicated by the multiple memcpy and memmove
  // std::vector<char> cigar_v(static_cast<unsigned long int>(qlen + dlen + 1));
  // char * cigar = cigar_v.data();
  char * cigar = static_cast<char*>(xmalloc(qlen + dlen + 1));

  char * cigarend {cigar + qlen + dlen + 1};

  auto operation = '\0';
  auto count = 0;
  *(--cigarend) = '\0';

  auto column = qlen;
  auto row = dlen;

  while ((column > 0) and (row > 0))
    {
      const auto cell = directions[qlen * (row - 1) + (column - 1)];

      ++alength;

      if ((operation == 'I') and ((cell & maskextleft) != 0))
        {
          --row;
          pushop('I', cigarend, operation, count);
        }
      else if ((operation == 'D') and ((cell & maskextup) != 0))
        {
          --column;
          pushop('D', cigarend, operation, count);
        }
      else if ((cell & maskleft) != 0)
        {
          --row;
          pushop('I', cigarend, operation, count);
        }
      else if ((cell & maskup) != 0)
        {
          --column;
          pushop('D', cigarend, operation, count);
        }
      else
        {
          if (nt_extract(qseq, column - 1) ==
              nt_extract(dseq, row - 1)) {
            ++matches;
          }
          --column;
          --row;
          pushop('M', cigarend, operation, count);
        }
    }

  while(column > 0)
    {
      ++alength;
      --column;
      pushop('D', cigarend, operation, count);
    }

  while(row > 0)
    {
      ++alength;
      --row;
      pushop('I', cigarend, operation, count);
    }

  finishop(cigarend, operation, count);

  /* move and reallocate cigar */

  auto cigaralloc = static_cast<std::size_t>(cigar + qlen + dlen - cigarend + 1);
  // note: std::memmove( void* dest, const void* src, std::size_t count );
  std::memmove(cigar, cigarend, cigaralloc);
  cigar = static_cast<char *>(xrealloc(cigar, cigaralloc));

  nwdiff = alength - matches;
  nwalignmentlength = alength;
  nwalignment = cigar;
}


/*

  Needleman/Wunsch/Sellers aligner

  finds a global alignment with minimum cost
  there should be positive costs/penalties for gaps and for mismatches
  matches should have zero cost (0)

  alignment priority when backtracking (from lower right corner):
  1. left/insert/e (gap in query sequence (qseq))
  2. align/diag/h (match/mismatch)
  3. up/delete/f (gap in database sequence (dseq))

  qseq: the reference/query/upper/vertical/from sequence
  dseq: the sample/database/lower/horizontal/to sequence

  matrix of qlen columns and dlen rows

  typical costs:
  match: 0
  mismatch: 3
  gapopen: 4
  gapextend: 3

  input

  dseq: pointer to start of database sequence
  dlen: database sequence length
  qseq: pointer to start of query sequence
  qlen: query sequence length
  score_matrix: 32x32 matrix of longs with scores for aligning two symbols
  gapopen: positive number indicating penalty for opening a gap of length zero
  gapextend: positive number indicating penalty for extending a gap

  output

  nwscore: the global alignment score
  nwdiff: number of non-identical nucleotides in one optimal global alignment
  nwalignmentlength: the length of one optimal alignment
  nwalignment: cigar string with one optimal alignment

*/

auto nw(char * dseq,
        const uint64_t dlen,
        char * qseq,
        const uint64_t qlen,
        const std::array<int64_t, n_cells * n_cells> & score_matrix,
        const uint64_t gapopen,
        const uint64_t gapextend,
        uint64_t & nwdiff,
        uint64_t & nwalignmentlength,
        char * & nwalignment,
        std::vector<unsigned char> & directions,
        std::vector<uint64_t> & hearray) -> void
{
  align(dseq, dlen, qseq, qlen, score_matrix,
          gapopen, gapextend, directions, hearray);

  backtrack(dseq, dlen, qseq, qlen, nwdiff,
            nwalignmentlength, nwalignment, directions);

  std::fill(directions.begin(), directions.end(), '\0');  // reset the alignment matrix
}
