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

#include "swarm.h"
#include <cstdint>


unsigned char * score_matrix_8 {nullptr};
unsigned short * score_matrix_16 {nullptr};
int64_t * score_matrix_63 {nullptr};


void score_matrix_read(struct Parameters const & p)
{
  static constexpr int n_cells {32};  // number of chars in sym_nt
  static constexpr long long int one_thousand {1000};
  long long int sc {0};
  long long int hi {-one_thousand};
  long long int lo {one_thousand};

  score_matrix_8 = new unsigned char[n_cells * n_cells];
  score_matrix_16 = new unsigned short[n_cells * n_cells];
  score_matrix_63 = new int64_t[n_cells * n_cells];

  for(auto a = 0; a < n_cells / 2; a++) {
    for(auto b = 0; b < n_cells / 2; b++) {
      sc = ((a == b) && (a > 0) && (b > 0)) ? 0 : p.penalty_mismatch;
      // sc = (a==b) ? matchscore : mismatchscore;
      if (sc < lo) {
        lo = sc;
      }
      if (sc > hi) {
        hi = sc;
      }
      score_matrix_63[(a * n_cells) + b] = sc;
    }
  }

  for(auto a = 0; a < n_cells; a++) {
    for(auto b = 0; b < n_cells; b++) {
      sc = score_matrix_63[(a * n_cells) + b];
      score_matrix_8[(a * n_cells) + b] = static_cast<unsigned char>(sc);
      score_matrix_16[(a * n_cells) + b] = static_cast<unsigned short>(sc);
    }
  }
}


void score_matrix_free()
{
  delete [] score_matrix_8;
  score_matrix_8 = nullptr;
  delete [] score_matrix_16;
  score_matrix_16 = nullptr;
  delete [] score_matrix_63;
  score_matrix_63 = nullptr;
}
