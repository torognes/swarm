/*
    SWARM

    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

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

unsigned char * score_matrix_8 = nullptr;
unsigned short * score_matrix_16 = nullptr;
int64_t * score_matrix_63 = nullptr;

void score_matrix_read();

#if 0

/* never used */

void score_matrix_dump()
{
  fprintf(logfile, "     ");
  for(auto i = 0; i < 16; i++)
    fprintf(logfile, "%2d", i);
  fprintf(logfile, "\n");
  fprintf(logfile, "     ");
  for(auto i = 0; i < 16; i++)
    fprintf(logfile, " %c", sym_nt[i]);
  fprintf(logfile, "\n");
  for(auto i = 0; i < 16; i++)
  {
    fprintf(logfile, "%2d %c ", i, sym_nt[i]);
    for(auto j = 0; j < 16; j++)
      {
        fprintf(logfile, "%2" PRId64, score_matrix_63[(i<<5) + j]);
      }
    fprintf(logfile, "\n");
  }
}

#endif

void score_matrix_read()
{
  constexpr int cells {32};
  constexpr long long int one_thousand {1000};
  constexpr unsigned int multiplier {5};
  long long int sc {0};
  long long int hi {-one_thousand};
  long long int lo {one_thousand};
  int64_t SCORELIMIT_8;
  int64_t SCORELIMIT_16;

  score_matrix_8 = static_cast<unsigned char *>(xmalloc(cells * cells * sizeof(char)));
  score_matrix_16 = static_cast<unsigned short *>(xmalloc(cells * cells * sizeof(short)));
  score_matrix_63 = static_cast<int64_t *>(xmalloc(cells * cells * sizeof(int64_t)));

  for(auto a = 0; a < cells / 2; a++) {
    for(auto b = 0; b < cells / 2; b++) {
      sc = ((a == b) && (a > 0) && (b > 0)) ? 0 : penalty_mismatch;
      // sc = (a==b) ? matchscore : mismatchscore;
      if (sc < lo) {
        lo = sc;
      }
      if (sc > hi) {
        hi = sc;
      }
      score_matrix_63[(a << multiplier) + b] = sc;
    }
  }


  SCORELIMIT_8  = UINT8_MAX  + 1 - hi;
  SCORELIMIT_16 = UINT16_MAX + 1 - hi;

  for(auto a = 0; a < cells; a++) {
    for(auto b = 0; b < cells; b++) {
      sc = score_matrix_63[(a << multiplier) + b];
      score_matrix_8[(a << multiplier) + b] = static_cast<unsigned char>(sc);
      score_matrix_16[(a << multiplier) + b] = static_cast<unsigned short>(sc);
    }
  }
}

void score_matrix_init()
{
  score_matrix_read();
  //  score_matrix_dump();
}

void score_matrix_free()
{
  xfree(score_matrix_8);
  score_matrix_8 = nullptr;
  xfree(score_matrix_16);
  score_matrix_16 = nullptr;
  xfree(score_matrix_63);
  score_matrix_63 = nullptr;
}
