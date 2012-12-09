/*
    SWARM

    Copyright (C) 2012 Torbjorn Rognes and Frederic Mahe

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

#define COMPDIFF
//#define SHOWALN

#ifdef COMPDIFF
BYTE diru[512*512];
BYTE dirl[512*512];
#endif

void nw(char * dseq,
	char * dend,
	char * qseq,
	char * qend,
	unsigned long * hearray,
	unsigned long * score_matrix,
	unsigned long gapopen,
	unsigned long gapextend,
	unsigned long * nwscore,
	unsigned long * nwdiff)
{
  unsigned long h, n, e, f;
  unsigned long *hep;
  char *qp, *dp;
  unsigned long * sp;

  dp = dseq;

#define MAXH (ULONG_MAX - 100)

#ifdef COMPDIFF
  memset(diru, 0, 512*512);
  memset(dirl, 0, 512*512);
#endif

  unsigned long qlen = qend - qseq;
  unsigned long dlen = dend - dseq;

  for(unsigned long i=0; i<2*qlen; i++)
  {
    hearray[2*i]   = gapopen + (i+1) * gapextend; // H 
    hearray[2*i+1] = gapopen + (i+1) * gapextend; // E
  }

  unsigned long i = 0;
  unsigned long j = 0;
  while (dp < dend)
    {
      hep = (unsigned long*) hearray;
      qp = qseq;
      sp = ((unsigned long*)(score_matrix)) + (*dp << 5);
      f = gapopen + (j+1) * gapextend;
      h = (j == 0) ? 0 : gapopen + j * gapextend;
      
      i = 0;
      while (qp < qend)
        {
          n = *hep;
          e = *(hep+1);
          h += sp[(unsigned)*qp];

          if (e < h)
            h = e;
          if (f < h)
            h = f;

#ifdef COMPDIFF
	  diru[qlen*j+i] = ((h == f) ? 1 : 0);
	  dirl[qlen*j+i] = ((h == e) ? 1 : 0);
#endif

          *hep = h;

          e += gapextend;
          f += gapextend;
          h += gapopen + gapextend;

          if (h < e)
            e = h;
          if (h < f)
            f = h;

          *(hep+1) = e;
          h = n;
          hep += 2;
          qp++;
	  i++;
        }

      dp++;
      j++;
    }

#ifdef COMPDIFF
  
#if 0

  for(i=1; i<=qlen; i++)
  {
    for(j=1; j<=dlen; j++)
    {
      if (diru[qlen*(j-1)+(i-1)])
      {
	if (dirl[qlen*(j-1)+(i-1)])
	  printf("+");
	else
	  printf("^");
      }
      else if (dirl[qlen*(j-1)+(i-1)])
      {
	printf("<");
      }
      else
      {
	printf("\\");
      }
    }
    printf("\n");
  }

#endif


#ifdef SHOWALN
  printf("Alignment: ");
#endif
  unsigned long diff = 0;

  i = qlen;
  j = dlen;
  
  while ((i>0) && (j>0))
  {
    if (diru[qlen*(j-1)+(i-1)])
    {
      diff++;
      i--;
#ifdef SHOWALN
      printf("D");
#endif
    }
    else if (dirl[qlen*(j-1)+(i-1)])
    {
      diff++;
      j--;
#ifdef SHOWALN
      printf("I");
#endif
    }
    else
    {
      if (qseq[i-1] == dseq[j-1])
      {
#ifdef SHOWALN
	printf("=");
#endif
      }
      else
      {
#ifdef SHOWALN
	printf("X");
#endif
	diff++;
      }
      i--;
      j--;
    }
  }
  
  while(i>0)
  {
    diff++;
    i--;
#ifdef SHOWALN
    printf("D");
#endif
  }
  
  while(j>0)
  {
    diff++;
    j--;
#ifdef SHOWALN
    printf("I");
#endif
  }

#ifdef SHOWALN
  printf("\n");
#endif

#endif

  unsigned long dist = hearray[2*qlen-2];

  //  printf("SLOW: dist: %lu  diff: %lu\n", dist, diff);
  
  * nwscore = dist;
  * nwdiff = diff;
}

