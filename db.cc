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

char map_nt[256] =
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  1, -1,  2, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  4,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1,  1, -1,  2, -1, -1, -1,  3, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1,  4,  4, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

char map_hex[256] =
  {
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, -1, -1, -1, -1, -1, -1,
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, 10, 11, 12, 13, 14, 15, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
    -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1
  };

unsigned long sequences = 0;
unsigned long nucleotides = 0;
unsigned long headerchars = 0;
int longest = 0;
int longestheader = 0;

seqinfo_t * seqindex = 0;
char * datap = 0;

#define MEMCHUNK 1048576
#define LINEALLOC LINE_MAX

void showseq(char * seq)
{
  char * p = seq;
  while (char c = *p++)
  {
    putchar(sym_nt[(unsigned int)c]);
  }
}


void db_read(const char * filename)
{
  /* allocate space */

  unsigned long dataalloc = MEMCHUNK;
  datap = (char *) xmalloc(dataalloc);
  unsigned long datalen = 0;

  longest = 0;
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;
  headerchars = 0;

  FILE * fp = NULL;
  if (filename)
    {
      fp = fopen(filename, "r");
      if (!fp)
	fatal("Error: Unable to open input data file (%s).", filename);
    }
  else
    fp = stdin;

  char line[LINEALLOC];
  line[0] = 0;
  fgets(line, LINEALLOC, fp);

  while(line[0])
    {
      /* read header */

      if (line[0] != '>')
        fatal("Illegal header line in fasta file.");
      
      long headerlen = strchrnul(line+1, '\n') - (line+1);

      headerchars += headerlen;

      if (headerlen > longestheader)
	longestheader = headerlen;


      /* store the header */

      while (datalen + headerlen + 1 > dataalloc)
      {
	dataalloc += MEMCHUNK;
	datap = (char *) xrealloc(datap, dataalloc);
      }
      memcpy(datap + datalen, line + 1, headerlen);
      *(datap + datalen + headerlen) = 0;
      datalen += headerlen + 1;


      /* get next line */

      line[0] = 0;
      fgets(line, LINEALLOC, fp);

      /* read sequence */

      unsigned long seqbegin = datalen;

      while (line[0] && (line[0] != '>'))
	{
	  char m;
	  char c;
	  char * p = line;
	  while((c = *p++))
	    if ((m = map_nt[(int)c]) >= 0)
	    {
	      while (datalen >= dataalloc)
	      {
		dataalloc += MEMCHUNK;
		datap = (char *) xrealloc(datap, dataalloc);
	      }

	      *(datap+datalen) = m;
	      datalen++;
	    }
	    else if (c != '\n')
	      fatal("Illegal character in sequence.");
	  line[0] = 0;
	  fgets(line, LINEALLOC, fp);
	}
      
      while (datalen >= dataalloc)
      {
	dataalloc += MEMCHUNK;
	datap = (char *) xrealloc(datap, dataalloc);
      }
      
      long length = datalen - seqbegin;

      nucleotides += length;

      if (length > longest)
	longest = length;

      *(datap+datalen) = 0;
      datalen++;

      sequences++;
    }

  fclose(fp);

  seqindex = (seqinfo_t *) xmalloc(sequences * sizeof(seqinfo_t));
  seqinfo_t * seqindex_p = seqindex;

  char * p = datap;
  for(unsigned long i=0; i<sequences; i++)
  {
    seqindex_p->header = p;
    seqindex_p->headerlen = strlen(p);
    seqindex_p->abundance = 0;
    sscanf(p, "%*[^_]_%lu\n", & seqindex_p->abundance);
    p += seqindex_p->headerlen + 1;

    /*
    fprintf(stdout, "Header: %s\n", seqindex_p->header);
    */

    seqindex_p->seq = p;
    seqindex_p->seqlen = strlen(p);
    p += seqindex_p->seqlen + 1;

    /*
    fprintf(stdout, "Sequence: ");
    showseq(seqindex_p->seq);
    fprintf(stdout, "\n");
    */

    seqindex_p++;
  }
}

unsigned long db_getsequencecount()
{
  return sequences;
}

unsigned long db_getnucleotidecount()
{
  return nucleotides;
}

unsigned long db_getlongestheader()
{
  return longestheader;
}

unsigned long db_getlongestsequence()
{
  return longest;
}

seqinfo_t * db_getseqinfo(unsigned long seqno)
{
  return seqindex+seqno;
}

char * db_getsequence(unsigned long seqno)
{
  return seqindex[seqno].seq;
}

void db_getsequenceandlength(unsigned long seqno,
			     char ** address,
			     long * length)
{
  *address = seqindex[seqno].seq;
  *length = (long)(seqindex[seqno].seqlen);
}

unsigned long db_getsequencelen(unsigned long seqno)
{
  return seqindex[seqno].seqlen;
}

char * db_getheader(unsigned long seqno)
{
  return seqindex[seqno].header;
}

unsigned long db_getheaderlen(unsigned long seqno)
{
  return seqindex[seqno].headerlen;
}

unsigned long db_getabundance(unsigned long seqno)
{
  return seqindex[seqno].abundance;
}

void db_putseq(long seqno)
{
  char * seq;
  long len;
  db_getsequenceandlength(seqno, & seq, & len);
  for(int i=0; i<len; i++)
    putchar(sym_nt[(int)(seq[i])]);
}


void db_showsequence(unsigned long seqno)
{
  unsigned long abundance, seqlen, hdrlen;
  char * seq, * hdr;

  abundance = db_getabundance(seqno);
  seqlen = db_getsequencelen(seqno);
  hdrlen = db_getheaderlen(seqno);
  seq = db_getsequence(seqno);
  hdr = db_getheader(seqno);

  printf(">%s\n", hdr);
  showseq(seq);
  printf("\n");
}

void db_showall()
{
  unsigned long seqcount = db_getsequencecount();

  for(unsigned long i = 0; i < seqcount; i++)
  {
    db_showsequence(i);
  }
}

void db_free()
{
  if (datap)
    free(datap);
  if (seqindex)
    free(seqindex);
}
