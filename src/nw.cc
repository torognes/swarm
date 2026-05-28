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

#include "nw.h"
#include "utils/cigar.h"
#include "utils/nt_codec.h"
#include <algorithm>  // std::min(), std::fill(), std::reverse()
#include <array>
#include <cassert>  // assert()
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // snprintf
#include <vector>


namespace {

constexpr unsigned char maskup      = 1;
constexpr unsigned char maskleft    = 2;
constexpr unsigned char maskextup   = 4;
constexpr unsigned char maskextleft = 8;


auto fill_matrix(char const * dseq,
                 const uint64_t dlen,
                 char const * qseq,
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
  //
  // value in a cell depends on three neighbors:
  // - cell on the 'left' (column - 1),
  // - cell above ('top' or 'up') (row - 1),
  // - 'diagonal' cell ('top' or 'up') (column - 1 and row - 1),
  static constexpr auto multiplier = 5U;

  assert(directions.size() >= qlen * dlen);
  assert(hearray.size() >= 2 * qlen);

  // hearray: array of alignments and insertions (initialized here,
  // then modified through a pointer; can't be const)
  for (auto column = 0UL; column < qlen; ++column) {
    hearray[2 * column]     = (1 * gapopen) + ((column + 1) * gapextend); // H (N)
    hearray[(2 * column) + 1] = (2 * gapopen) + ((column + 2) * gapextend); // E
  }

  uint64_t previous_diagonal {0};
  uint64_t left {0};

  for (auto row = 0UL; row < dlen; ++row) {
      auto he_index = 0UL;
      auto top = (2 * gapopen) + ((row + 2) * gapextend);
      uint64_t diagonal = (row == 0) ? 0 : (gapopen + (row * gapextend));

      for (auto column = 0UL; column < qlen; ++column) {
          const auto index = (qlen * row) + column;

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


auto backtrack(char const * dseq,
               const uint64_t dlen,
               char const * qseq,
               const uint64_t qlen,
               uint64_t & nwdiff,
               std::vector<unsigned char> const & directions,
               std::vector<char> & raw_alignment) -> void
{
  /* backtrack: count differences and save alignment in cigar string */

  uint64_t alength {0};  // refactoring: eliminate variable
  uint64_t matches {0};

  auto operation = '\0';

  auto column = qlen;
  auto row = dlen;

  while ((column > 0) and (row > 0))
    {
      const auto cell = directions[(qlen * (row - 1)) + (column - 1)];

      ++alength;

      if ((operation == 'I') and ((cell & maskextleft) != 0))
        {
          --row;
          raw_alignment.emplace_back('I');
          operation = 'I';
        }
      else if ((operation == 'D') and ((cell & maskextup) != 0))
        {
          --column;
          raw_alignment.emplace_back('D');
          operation = 'D';
        }
      else if ((cell & maskleft) != 0)
        {
          --row;
          raw_alignment.emplace_back('I');
          operation = 'I';
        }
      else if ((cell & maskup) != 0)
        {
          --column;
          raw_alignment.emplace_back('D');
          operation = 'D';
        }
      else
        {
          if (nt_extract(qseq, column - 1) == nt_extract(dseq, row - 1)) {
            ++matches;
          }
          --column;
          --row;
          raw_alignment.emplace_back('M');
          operation = 'M';
        }
    }

  while (column > 0) {
      ++alength;
      --column;
      raw_alignment.emplace_back('D');
    }

  while (row > 0) {
      ++alength;
      --row;
      raw_alignment.emplace_back('I');
    }

  nwdiff = alength - matches;
  assert(raw_alignment.size() == alength);
}

}  // unnamed namespace


NwAligner::NwAligner(uint64_t const longest_sequence,
                     int64_t const penalty_mismatch,
                     uint64_t const gapopen,
                     uint64_t const gapextend)
  : directions_(longest_sequence * longest_sequence),
    hearray_(2 * longest_sequence),
    score_matrix_(create_score_matrix<int64_t>(penalty_mismatch)),
    gapopen_(gapopen),
    gapextend_(gapextend)
{
  raw_alignment_.reserve(2 * longest_sequence);
  cigar_string_.reserve(2 * longest_sequence);
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

  input (per align() call)

  dseq: pointer to start of database sequence
  dlen: database sequence length
  qseq: pointer to start of query sequence
  qlen: query sequence length

  scoring configuration (set once at NwAligner construction)

  score_matrix: 32x32 matrix of longs with scores for aligning two symbols
                (derived from penalty_mismatch)
  gapopen: positive number indicating penalty for opening a gap of length zero
  gapextend: positive number indicating penalty for extending a gap

  output (carried in NwAligner::Result)

  nwscore: the global alignment score
  nwdiff: number of non-identical nucleotides in one optimal global alignment
          (Result::differences)
  nwalignmentlength: the length of one optimal alignment (Result::length)
  nwalignment: cigar string with one optimal alignment (Result::cigar_string)

*/

auto NwAligner::align(char const * dseq, uint64_t const dlen,
                      char const * qseq, uint64_t const qlen) -> NwAligner::Result
{
  static constexpr auto one_hundred = 100.0;

  raw_alignment_.clear();
  cigar_string_.clear();

  fill_matrix(dseq, dlen, qseq, qlen, score_matrix_,
              gapopen_, gapextend_, directions_, hearray_);

  uint64_t nwdiff {0};
  backtrack(dseq, dlen, qseq, qlen, nwdiff, directions_, raw_alignment_);

  std::fill(directions_.begin(), directions_.end(), '\0');  // reset the alignment matrix

  // backtracking produces a reversed alignment (starting from the end)
  std::reverse(raw_alignment_.begin(), raw_alignment_.end());
  compress_alignment_to_cigar(raw_alignment_, cigar_string_);

  // loosing precision when converting raw_alignment_.size() and nwdiff
  // to double is not an issue, no need to add assertions
  auto const length = raw_alignment_.size();
  auto const percent_id =
    one_hundred * static_cast<double>(length - nwdiff) / static_cast<double>(length);

  return Result{cigar_string_, nwdiff, length, percent_id};
}
