/*
  SWARM

  Copyright (C) 2012-2014 Torbjorn Rognes and Frederic Mahe

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

//#define HASH hash_fnv_1a_64
#define HASH hash_cityhash64

#define MEMCHUNK 1048576
#define LINEALLOC LINE_MAX

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
qgramvector_t * qgrams = 0;

void showseq(char * seq)
{
  char * p = seq;
  while (char c = *p++)
    {
      putchar(sym_nt[(unsigned int)c]);
    }
}

void fprint_id(FILE * stream, unsigned long x)
{
  fprintf(stream, "%.*s", seqindex[x].headeridlen, seqindex[x].header);
}

void fprint_id_noabundance(FILE * stream, unsigned long x)
{
  seqinfo_t * sp = seqindex + x;
  char * h = sp->header;
  int hdrlen = sp->headerlen;

  /* print start of header */
  fprintf(stream, "%.*s", sp->abundance_start, h);

  /* print semicolon if the abundance is not at the end (usearch style) */
  if ((sp->abundance_start > 0) && (sp->abundance_end < hdrlen))
    fprintf(stream, ";");

  /* print remaining part */
  fprintf(stream, "%.*s", hdrlen - sp->abundance_end, h + sp->abundance_end);
}

int db_compare_abundance(const void * a, const void * b)
{
  seqinfo_t * x = (seqinfo_t *) a;
  seqinfo_t * y = (seqinfo_t *) b;
  
  if (x->abundance > y->abundance)
    return -1;
  else if (x->abundance < y->abundance)
    return +1;
#if 1
  else if (x < y)
    return -1;
  else if (x > y)
    return +1;
  else
    return 0;
#else
  else
    return strcmp(x->seq, y->seq);
#endif
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

  /* get file size */

  if (fseek(fp, 0, SEEK_END))
    fatal("Error: Unable to seek in database file (%s)", filename);

  long filesize = ftell(fp);
  
  rewind(fp);

  char line[LINEALLOC];
  line[0] = 0;
  fgets(line, LINEALLOC, fp);

  long lineno = 1;

  progress_init("Reading database: ", filesize);
  while(line[0])
    {
      /* read header */

      if (line[0] != '>')
        fatal("Illegal header line in fasta file.");
      
      long headerlen = xstrchrnul(line+1, '\n') - (line+1);

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
      lineno++;


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
              {
                char msg[100];
                snprintf(msg, 100, "Illegal character '%c' in sequence on line %ld", c, lineno);
                fatal(msg);
              }
          line[0] = 0;
          fgets(line, LINEALLOC, fp);
          lineno++;
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
      
      progress_update(ftell(fp));
    }
  progress_done();

  fclose(fp);

  /* set up hash to check for unique headers */

  unsigned long hdrhashsize = 2 * sequences;

  seqinfo_t * * hdrhashtable = 
    (seqinfo_t **) xmalloc(hdrhashsize * sizeof(seqinfo_t *));
  memset(hdrhashtable, 0, hdrhashsize * sizeof(seqinfo_t *));

  unsigned long duplicatedidentifiers = 0;
  
  /* create indices */

  seqindex = (seqinfo_t *) xmalloc(sequences * sizeof(seqinfo_t));
  seqinfo_t * seqindex_p = seqindex;

  regex_t db_regexp;
  regmatch_t pmatch[4];

  if (usearch_abundance)
    if (regcomp(&db_regexp, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
      fatal("Regular expression compilation failed");

  long lastabundance = LONG_MAX;

  int presorted = 1;

  char * p = datap;
  progress_init("Indexing database:", sequences);
  for(unsigned long i=0; i<sequences; i++)
    {
      seqindex_p->header = p;
      seqindex_p->headerlen = strlen(seqindex_p->header);
      p += seqindex_p->headerlen + 1;

      seqindex_p->headeridlen = xstrchrnul(seqindex_p->header, ' ') 
        - seqindex_p->header;

      /* get amplicon abundance */
      seqindex_p->abundance = 1;
      if (usearch_abundance)
        {
          if (!regexec(&db_regexp, seqindex_p->header, 4, pmatch, 0))
            {
              seqindex_p->abundance = atol(seqindex_p->header + pmatch[2].rm_so);
              seqindex_p->abundance_start = pmatch[0].rm_so;
              seqindex_p->abundance_end = pmatch[0].rm_eo;
            }
          else
            {
              seqindex_p->abundance_start = 0;
              seqindex_p->abundance_end = 0;
            }
        }
      else
        {
          int start, end;
          unsigned long abundance;
          if (sscanf(seqindex_p->header, "%*[^ _]%n_%lu%n",
                     & start,
                     & abundance,
                     & end) == 1)
            {
              seqindex_p->abundance = abundance;
              seqindex_p->abundance_start = start;
              seqindex_p->abundance_end = end;
            }
          else
            {
              seqindex_p->abundance_start = 0;
              seqindex_p->abundance_end = 0;
            }
        }

      if (seqindex_p->abundance > lastabundance)
        presorted = 0;

      lastabundance = seqindex_p->abundance;

      if (seqindex_p->abundance < 1)
        fatal("abundance annotation zero");

      /* check hash, fatal error if found, otherwize insert new */
      unsigned long hdrhash = HASH((unsigned char*)seqindex_p->header, seqindex_p->headeridlen);
      seqindex_p->hdrhash = hdrhash;
      unsigned long hashindex = hdrhash % hdrhashsize;

      seqinfo_t * found;
    
      while ((found = hdrhashtable[hashindex]))
        {
          if ((found->hdrhash == hdrhash) &&
              (found->headeridlen == seqindex_p->headeridlen) &&
              (strncmp(found->header, seqindex_p->header, found->headeridlen) == 0))
            break;
          hashindex = (hashindex + 1) % hdrhashsize;
        }

      if (found)
        {
          duplicatedidentifiers++;
          fprintf(stderr, "Duplicated sequence identifier: %s\n", seqindex_p->header);
        }

      hdrhashtable[hashindex] = seqindex_p;
    
      seqindex_p->seq = p;
      seqindex_p->seqlen = strlen(p);
      p += seqindex_p->seqlen + 1;

      seqindex_p++;
      progress_update(i);
    }
  progress_done();

#if 1
  if (!presorted)
#endif
    {
      progress_init("Abundance sorting:", 1);
      qsort(seqindex, sequences, sizeof(seqinfo_t), db_compare_abundance);
      progress_done();
    }

  if (usearch_abundance)
    regfree(&db_regexp);

  free(hdrhashtable);

  if (duplicatedidentifiers)
    exit(1);
}

void db_qgrams_init()
{
  qgrams = (qgramvector_t *) xmalloc(sequences * sizeof(qgramvector_t));

  seqinfo_t * seqindex_p = seqindex;
  progress_init("Find qgram vects: ", sequences);
  for(int i=0; i<sequences; i++)
    {
      /* find qgrams */
      findqgrams((unsigned char*) seqindex_p->seq,
                 seqindex_p->seqlen,
                 qgrams[i]);
      seqindex_p++;
      progress_update(i);
    }
  progress_done();
}

void db_qgrams_done()
{
  free(qgrams);
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

void db_free()
{
  if (datap)
    free(datap);
  if (seqindex)
    free(seqindex);
}
