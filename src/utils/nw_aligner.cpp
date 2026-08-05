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

#include "nw_aligner.hpp"
#include "../db.hpp"  // struct Sequence
#include "cigar.hpp"
#include "score_matrix.hpp"  // n_cells, create_score_matrix
#include "span.hpp"  // Span, make_span
#include "view.hpp"  // View, make_view
#include <algorithm>  // std::min(), std::reverse()
#include <array>
#include <cassert>  // assert()
#include <cstdint>  // int64_t, uint64_t
#include <vector>


namespace {

constexpr unsigned char maskup      = 1;
constexpr unsigned char maskleft    = 2;
constexpr unsigned char maskextup   = 4;
constexpr unsigned char maskextleft = 8;


auto fill_matrix(Sequence const & dseq,
                 Sequence const & qseq,
                 const std::array<int64_t, n_cells * n_cells> & score_matrix,
                 const uint64_t gapopen,
                 const uint64_t gapextend,
                 Span<unsigned char> const directions,
                 Span<NwAligner::HECell> const hearray) -> void
{
  // Sequence::length is the nucleotide count (not encoded.size(), which
  // is the packed-byte count); nucleotide_at() and the inner loops below
  // both work in nucleotide units.
  auto const dlen = static_cast<uint64_t>(dseq.length);
  auto const qlen = static_cast<uint64_t>(qseq.length);

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
  auto const new_gap = gapopen + gapextend;

  assert(directions.size() >= qlen * dlen);
  assert(hearray.size() >= qlen);

  // hearray: per-column carry of (H, E) updated as rows advance.
  // This loop seeds the row-0 boundary; subsequent rows overwrite in
  // place.
  for (auto column = 0UL; column < qlen; ++column) {
    hearray[column].h_score = gapopen + ((column + 1) * gapextend);
    hearray[column].e_score = (2 * gapopen) + ((column + 2) * gapextend);
  }

  for (auto row = 0UL; row < dlen; ++row) {
      auto top = (2 * gapopen) + ((row + 2) * gapextend);
      uint64_t diagonal = (row == 0) ? 0 : (gapopen + (row * gapextend));
      auto const row_offset = (nucleotide_at(dseq, row) + 1U) << multiplier;

      // this row of the direction matrix, so the inner loop indexes a row
      // rather than recomputing (qlen * row) + column into the whole
      // matrix, and its writes are bounds-checked against the row
      auto const row_directions = directions.subspan(qlen * row, qlen);

      for (auto column = 0UL; column < qlen; ++column) {
          auto const previous_diagonal = hearray[column].h_score;
          auto left                    = hearray[column].e_score;
          unsigned char flags          = '\0';

          diagonal += static_cast<uint64_t>(
              score_matrix[row_offset + nucleotide_at(qseq, column) + 1U]);

          flags |= (top < diagonal) ? maskup : 0U;
          diagonal = std::min({diagonal, top, left});
          flags |= (left == diagonal) ? maskleft : 0U;

          hearray[column].h_score = diagonal;

          diagonal += new_gap;
          left += gapextend;
          top  += gapextend;

          flags |= (top  < diagonal) ? maskextup   : 0U;
          flags |= (left < diagonal) ? maskextleft : 0U;
          top  = std::min(diagonal, top);
          left = std::min(diagonal, left);

          row_directions[column]  = flags;
          hearray[column].e_score = left;
          diagonal                = previous_diagonal;
        }
    }
}


// directions is indexed flat here rather than a row at a time: the walk
// moves up and left, so it leaves a row as often as it stays in one, and a
// per-step subview would cost more than the multiplication it saves.
auto backtrack(Sequence const & dseq,
               Sequence const & qseq,
               View<unsigned char> const directions,
               std::vector<char> & raw_alignment) -> uint64_t
{
  /* backtrack: count differences and save alignment in cigar string */

  auto const dlen = static_cast<uint64_t>(dseq.length);
  auto const qlen = static_cast<uint64_t>(qseq.length);

  uint64_t matches {0};

  auto operation = '\0';

  auto column = qlen;
  auto row = dlen;

  while ((column > 0) and (row > 0))
    {
      const auto cell = directions[(qlen * (row - 1)) + (column - 1)];

      const bool cell_extends_left = (cell & maskextleft) != 0;
      const bool cell_extends_up   = (cell & maskextup)   != 0;
      const bool cell_opens_left   = (cell & maskleft)    != 0;
      const bool cell_opens_up     = (cell & maskup)      != 0;

      // Priority: extending an ongoing gap (sticky) beats opening a
      // new one; among fresh openings, insert beats delete; otherwise
      // take the diagonal (match/mismatch).
      const bool extend_insert = (operation == 'I') and cell_extends_left;
      const bool extend_delete = (operation == 'D') and cell_extends_up;
      const bool is_insert     = extend_insert or (not extend_delete and cell_opens_left);
      const bool is_delete     = extend_delete or cell_opens_up;

      if (is_insert)
        {
          --row;
          raw_alignment.emplace_back('I');
          operation = 'I';
        }
      else if (is_delete)
        {
          --column;
          raw_alignment.emplace_back('D');
          operation = 'D';
        }
      else
        {
          if (nucleotide_at(qseq, column - 1) == nucleotide_at(dseq, row - 1)) {
            ++matches;
          }
          --column;
          --row;
          raw_alignment.emplace_back('M');
          operation = 'M';
        }
    }

  // emit the remaining boundary: any unconsumed qseq columns become
  // deletions, any unconsumed dseq rows become insertions
  raw_alignment.insert(raw_alignment.end(), column, 'D');
  raw_alignment.insert(raw_alignment.end(), row,    'I');

  return raw_alignment.size() - matches;
}

}  // unnamed namespace


NwAligner::NwAligner(uint64_t const longest_sequence,
                     int64_t const penalty_mismatch,
                     uint64_t const gapopen,
                     uint64_t const gapextend)
  : directions_(longest_sequence * longest_sequence),
    hearray_(longest_sequence),
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

  dseq: database sequence (packed bytes + nucleotide length)
  qseq: query sequence    (packed bytes + nucleotide length)

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

auto NwAligner::align(Sequence const & dseq, Sequence const & qseq) -> NwAligner::Result
{
  static constexpr auto one_hundred = 100.0;

  raw_alignment_.clear();
  cigar_string_.clear();

  // fill_matrix() writes every directions[i] for i in [0, dlen*qlen),
  // so the buffer's content on entry doesn't matter.
  fill_matrix(dseq, qseq, score_matrix_,
              gapopen_, gapextend_,
              make_span(directions_), make_span(hearray_));

  auto const nwdiff = backtrack(dseq, qseq, make_view(directions_), raw_alignment_);

  // backtracking produces a reversed alignment (starting from the end)
  std::reverse(raw_alignment_.begin(), raw_alignment_.end());
  compress_alignment_to_cigar(make_view(raw_alignment_), cigar_string_);

  // loosing precision when converting raw_alignment_.size() and nwdiff
  // to double is not an issue, no need to add assertions
  auto const length = raw_alignment_.size();
  // an empty alignment would make percent_id a 0/0 NaN below; swarm
  // rejects empty input sequences upstream (db.cpp, "Empty sequence
  // found"), so the alignment of two non-empty sequences is non-empty.
  assert(length != 0);
  auto const percent_id =
    one_hundred * static_cast<double>(length - nwdiff) / static_cast<double>(length);

  return Result{make_view(cigar_string_), nwdiff, length, percent_id};
}
