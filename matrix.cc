/*
    SWARM

    Copyright (C) 2012-2013 Torbjorn Rognes and Frederic Mahe

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

long SCORELIMIT_7 = 0;
long SCORELIMIT_8;
long SCORELIMIT_16;
long SCORELIMIT_32;
long SCORELIMIT_63;
char BIAS;

unsigned char * score_matrix_8 = NULL;
unsigned short * score_matrix_16 = NULL;
long * score_matrix_63 = NULL;

void score_matrix_dump()
{
  fprintf(stderr, "     ");
  for(int i=0; i<16; i++)
    fprintf(stderr, "%2d", i);
  fprintf(stderr, "\n");
  fprintf(stderr, "     ");
  for(int i=0; i<16; i++)
    fprintf(stderr, " %c", sym_nt[i]);
  fprintf(stderr, "\n");
  for(int i=0; i<16; i++)
  {
    fprintf(stderr, "%2d %c ", i, sym_nt[i]);
    for(int j=0; j<16; j++)
      {
	fprintf(stderr, "%2ld", score_matrix_63[(i<<5) + j]);
      }
    fprintf(stderr, "\n");
  }
}

void score_matrix_read()
{
  int a, b;
  long sc, lo, hi; 
  
  score_matrix_8 = (unsigned char *) xmalloc(32*32*sizeof(char));
  score_matrix_16 = (unsigned short *) xmalloc(32*32*sizeof(short));
  score_matrix_63 = (long *) xmalloc(32*32*sizeof(long));
  
  hi = -1000;
  lo = 1000;

  for(a=0;a<16;a++)
    for(b=0;b<16;b++)
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
  
  for(a=0;a<32;a++)
    for(b=0;b<32;b++)
    {
      sc = score_matrix_63[(a<<5) + b];
      score_matrix_8[(a<<5) + b] = (unsigned char) sc;
      score_matrix_16[(a<<5) + b] = (unsigned short) sc;
    }
}

void score_matrix_init()
{
  score_matrix_read();
  //  score_matrix_dump();
}

void score_matrix_free()
{
  free(score_matrix_8);
  score_matrix_8 = NULL;
  free(score_matrix_16);
  score_matrix_16 = NULL;
  free(score_matrix_63);
  score_matrix_63 = NULL;
}
