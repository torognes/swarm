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

void nw(char * dseq,
	char * dend,
	char * qseq,
	char * qend,
	long * score_matrix,
	unsigned long gapopen,
	unsigned long gapextend,
	unsigned long * nwscore,
	unsigned long * nwdiff,
	unsigned long * nwalignmentlength,
	char ** nwalignment)
{
  unsigned long h, n, e, f;
  unsigned long *hep;
  char *qp, *dp;
  unsigned long * sp;

  unsigned long qlen = qend - qseq;
  unsigned long dlen = dend - dseq;

  BYTE * diru = (BYTE*) xmalloc(qlen*dlen);
  BYTE * dirl = (BYTE*) xmalloc(qlen*dlen);

  memset(diru, 0, qlen*dlen);
  memset(dirl, 0, qlen*dlen);

  unsigned long * hearray = 
    (unsigned long *) xmalloc(2 * qlen * sizeof(unsigned long));

  for(unsigned long i=0; i<qlen; i++)
  {
    hearray[2*i]   = gapopen + (i+1) * gapextend; // H 
    hearray[2*i+1] = gapopen + (i+1) * gapextend; // E
  }

  unsigned long i = 0;
  unsigned long j = 0;

  dp = dseq;

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

	  diru[qlen*j+i] = ((h == f) ? 1 : 0);
	  dirl[qlen*j+i] = ((h == e) ? 1 : 0);

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

  unsigned long dist = hearray[2*qlen-2];

  free(hearray);

  unsigned long diff = 0;

  /* backtrack: count differences and save alignment */

  i = qlen;
  j = dlen;

  char * alignmentstring = (char *) xmalloc(qlen + dlen + 1);
  char * p = alignmentstring;
  
  while ((i>0) && (j>0))
  {
    if (diru[qlen*(j-1)+(i-1)])
    {
      diff++;
      i--;
      *p++ = 'D';
    }
    else if (dirl[qlen*(j-1)+(i-1)])
    {
      diff++;
      j--;
      *p++ = 'I';
    }
    else
      {
	if (qseq[i-1] != dseq[j-1])
	  diff++;
	i--;
	j--;
	*p++ = 'M';
      }
  }
  
  while(i>0)
  {
    diff++;
    i--;
    *p++ = 'D';
  }
  
  while(j>0)
  {
    diff++;
    j--;
    *p++ = 'I';
  }

  *p = 0;

  free(diru);
  free(dirl);

  unsigned long alignmentlength = p - alignmentstring;

  /* reverse and compress alignment string */

  char * compressedalignment = (char *) xmalloc(alignmentlength + 1);
  char * q = compressedalignment;

  char x = 0;
  unsigned long no = 0;
  while (p-- > alignmentstring)
    {
      char c = *p;
      if (c == x)
	no++;
      else
	{
	  if (x)
	    {
	      if (no>1)
		q += sprintf(q, "%lu%c", no, x);
	      else
		*q++ = x;
	    }
	  x = c;
	  no = 1;
	}
    }
  
  if (x)
    {
      if (no>1)
	q += sprintf(q, "%lu%c", no, x);
      else
	*q++ = x;
    }
  
  *q = 0;
  
  free(alignmentstring);
  
  unsigned long compressedalignmentlength = q - compressedalignment;
  compressedalignment = (char*) xrealloc(compressedalignment, compressedalignmentlength + 1);

  * nwscore = dist;
  * nwdiff = diff;
  * nwalignmentlength = alignmentlength;
  * nwalignment = compressedalignment;
}
