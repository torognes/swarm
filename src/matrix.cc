/*
    SWARM

    Copyright (C) 2012-2020 Torbjorn Rognes and Frederic Mahe

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

int64_t SCORELIMIT_7 = 0;
int64_t SCORELIMIT_8;
int64_t SCORELIMIT_16;
int64_t SCORELIMIT_32;
int64_t SCORELIMIT_63;
char BIAS;

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
  int64_t sc, lo, hi;

  score_matrix_8 = static_cast<unsigned char*>(xmalloc(32*32*sizeof(char)));
  score_matrix_16 = static_cast<unsigned short*>(xmalloc(32*32*sizeof(short)));
  score_matrix_63 = static_cast<int64_t *>(xmalloc(32*32*sizeof(int64_t)));

  hi = -1000;
  lo = 1000;

  for(auto a = 0; a < 16; a++)
    for(auto b = 0; b < 16; b++)
    {
      sc = ((a==b)&&(a>0)&&(b>0)) ? 0 : penalty_mismatch;
      // sc = (a==b) ? matchscore : mismatchscore;
      if (sc < lo)
        lo = sc;
      if (sc > hi)
        hi = sc;
      score_matrix_63[(a<<5) + b] = sc;
    }


  SCORELIMIT_8  = 256 - hi;
  SCORELIMIT_16 = 65536 - hi;

  for(auto a = 0; a < 32; a++)
    for(auto b = 0; b < 32; b++)
    {
      sc = score_matrix_63[(a<<5) + b];
      score_matrix_8[(a<<5) + b] = static_cast<unsigned char>(sc);
      score_matrix_16[(a<<5) + b] = static_cast<unsigned short>(sc);
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
