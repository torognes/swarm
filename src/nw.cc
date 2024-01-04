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

auto pushop(const char newop, char * & cigarendp, char & op, int & count) -> void
{
  static constexpr auto buffer_length = 25U;

  if (newop == op) {
    ++count;
  }
  else
  {
    --cigarendp;
    *cigarendp = op;  // write char op at position end - 1
    if (count > 1)
    {
      std::array<char, buffer_length> buf {{}};
      const auto len = std::snprintf(buf.data(), buffer_length, "%d", count);
      assert(len >= 0);
      cigarendp -= len;
      std::memcpy(cigarendp, buf.data(), static_cast<std::size_t>(len));
    }
    op = newop;
    count = 1;
  }
}

auto finishop(char * & cigarendp, char & op, int & count) -> void
{
  static constexpr auto buffer_length = 25U;

  if ((op != '\0') and (count != 0))
  {
    --cigarendp;
    *cigarendp = op;  // write char op at position end - 1
    if (count > 1)
    {
      std::array<char, buffer_length> buf {{}};
      const auto len = std::snprintf(buf.data(), buffer_length, "%d", count);
      assert(len >= 0);
      cigarendp -= len;
      std::memcpy(cigarendp, buf.data(), static_cast<std::size_t>(len));
    }
    op = '\0';
    count = 0;
  }
}

constexpr unsigned char maskup      = 1;
constexpr unsigned char maskleft    = 2;
constexpr unsigned char maskextup   = 4;
constexpr unsigned char maskextleft = 8;

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
        uint64_t & nwscore,
        uint64_t & nwdiff,
        uint64_t & nwalignmentlength,
        char ** nwalignment,
        std::vector<unsigned char> & dir,
        std::vector<uint64_t> & hearray) -> void
{
  assert(dir.size() >= qlen * dlen);
  assert(hearray.size() >= 2 * qlen);
  static constexpr auto multiplier = 5U;

  uint64_t n {0};
  uint64_t e {0};

  for(auto i = 0UL; i < qlen; i++)
    {
      hearray[2 * i]     = 1 * gapopen + (i + 1) * gapextend; // H (N)
      hearray[2 * i + 1] = 2 * gapopen + (i + 2) * gapextend; // E
    }

  for(auto j = 0UL; j < dlen; j++)
    {
      auto * hep = hearray.data();
      auto f = 2 * gapopen + (j + 2) * gapextend;
      uint64_t h = (j == 0) ? 0 : (gapopen + j * gapextend);

      for(auto i = 0UL; i < qlen; i++)
        {
          const auto index = qlen * j + i;

          n = *hep;
          e = *(hep + 1);
          h += static_cast<uint64_t>(score_matrix
                                     [((nt_extract(dseq, j) + 1U) << multiplier)
                                      + (nt_extract(qseq, i) + 1)]);

          dir[index] |= (f < h ? maskup : 0U);
          h = std::min(h, f);
          h = std::min(h, e);
          dir[index] |= (e == h ? maskleft : 0U);

          *hep = h;

          h += gapopen + gapextend;
          e += gapextend;
          f += gapextend;

          dir[index] |= (f < h ? maskextup : 0U);
          dir[index] |= (e < h ? maskextleft : 0U);
          f = std::min(h, f);
          e = std::min(h, e);

          *(hep + 1) = e;
          h = n;
          hep += 2;
        }
    }

  const auto dist = hearray[2 * qlen - 2];

  /* backtrack: count differences and save alignment in cigar string */
  uint64_t alength {0};
  uint64_t matches {0};

  // refactoring: complicated by the multiple memcpy and memmove
  // std::vector<char> cigar_v(static_cast<unsigned long int>(qlen + dlen + 1));
  // char * cigar = cigar_v.data();
  char * cigar = static_cast<char*>(xmalloc(qlen + dlen + 1));

  char * cigarend {cigar + qlen + dlen + 1};

  auto op = '\0';
  auto count = 0;
  *(--cigarend) = '\0';

  auto i = qlen;
  auto j = dlen;

  while ((i > 0) and (j > 0))
    {
      const auto d = dir[qlen * (j - 1) + (i - 1)];  // refactoring: how to rename?

      ++alength;

      if ((op == 'I') and ((d & maskextleft) != 0))
        {
          --j;
          pushop('I', cigarend, op, count);
        }
      else if ((op == 'D') and ((d & maskextup) != 0))
        {
          --i;
          pushop('D', cigarend, op, count);
        }
      else if ((d & maskleft) != 0)
        {
          --j;
          pushop('I', cigarend, op, count);
        }
      else if ((d & maskup) != 0)
        {
          --i;
          pushop('D', cigarend, op, count);
        }
      else
        {
          if (nt_extract(qseq, i - 1) ==
              nt_extract(dseq, j - 1)) {
            ++matches;
          }
          --i;
          --j;
          pushop('M', cigarend, op, count);
        }
    }

  while(i > 0)
    {
      ++alength;
      --i;
      pushop('D', cigarend, op, count);
    }

  while(j > 0)
    {
      ++alength;
      --j;
      pushop('I', cigarend, op, count);
    }

  finishop(cigarend, op, count);

  /* move and reallocate cigar */

  auto cigaralloc = static_cast<std::size_t>(cigar + qlen + dlen - cigarend + 1);
  // note: std::memmove( void* dest, const void* src, std::size_t count );
  std::memmove(cigar, cigarend, cigaralloc);
  cigar = static_cast<char *>(xrealloc(cigar, cigaralloc));

  nwscore = dist;
  nwdiff = alength - matches;
  nwalignmentlength = alength;
  * nwalignment = cigar;
  std::fill(dir.begin(), dir.end(), '\0');  // reset the alignment matrix
}
