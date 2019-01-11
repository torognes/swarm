/*
    SWARM

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

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
  seqinfo_t * sp = seqindex + x;
  char * h = sp->header;
  int hdrlen = sp->headerlen;

  if (opt_append_abundance && (sp->abundance_start == sp->abundance_end))
    if (opt_usearch_abundance)
      fprintf(stream, "%.*s;size=%lu;", hdrlen, h, sp->abundance);
    else
      fprintf(stream, "%.*s_%lu", hdrlen, h, sp->abundance);
  else
    fprintf(stream, "%.*s", hdrlen, h);
}

void fprint_id_noabundance(FILE * stream, unsigned long x)
{
  seqinfo_t * sp = seqindex + x;
  char * h = sp->header;
  int hdrlen = sp->headerlen;

  if (sp->abundance_start < sp->abundance_end)
    {
      /* print start of header */
      fprintf(stream, "%.*s", sp->abundance_start, h);
      
      if (opt_usearch_abundance)
        {
          /* print semicolon if the abundance is not at either end */
          if ((sp->abundance_start > 0) && (sp->abundance_end < hdrlen))
            fprintf(stream, ";");
          
          /* print remaining part */
          fprintf(stream, "%.*s", hdrlen - sp->abundance_end, h + sp->abundance_end);
        }
    }
  else
    fprintf(stream, "%.*s", hdrlen, h);
}

void fprint_id_with_new_abundance(FILE * stream,
                                  unsigned long seqno,
                                  unsigned long abundance)
{
  seqinfo_t * sp = seqindex + seqno;

  if (opt_usearch_abundance)
    fprintf(stream,
            "%.*s%ssize=%lu;%.*s",
            sp->abundance_start,
            sp->header,
            sp->abundance_start > 0 ? ";" : "",
            abundance,
            sp->headerlen - sp->abundance_end,
            sp->header + sp->abundance_end);
  else
    fprintf(stream,
            "%.*s_%lu",
            sp->abundance_start,
            sp->header,
            abundance);
}

int db_compare_abundance(const void * a, const void * b)
{
  seqinfo_t * x = (seqinfo_t *) a;
  seqinfo_t * y = (seqinfo_t *) b;
  
  if (x->abundance > y->abundance)
    return -1;
  else if (x->abundance < y->abundance)
    return +1;
  else 
    return strcmp(x->header, y->header);
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

  struct stat fs;

  if (fstat(fileno(fp), & fs))
    fatal("Unable to fstat on input file (%s)", filename);
  bool is_regular = S_ISREG(fs.st_mode);
  long filesize = is_regular ? fs.st_size : 0;

  if (! is_regular)
    fprintf(logfile, "Waiting for data... (Hit Ctrl-C and run swarm -h if you meant to read data from a file.)\n");

  char line[LINEALLOC];
  line[0] = 0;
  if (!fgets(line, LINEALLOC, fp))
    line[0] = 0;

  unsigned int lineno = 1;

  progress_init("Reading database: ", filesize);
  while(line[0])
    {
      /* read header */
      /* the header ends at a space character, a newline or a nul character */

      if (line[0] != '>')
        fatal("Illegal header line in fasta file.");
      
      long headerlen = 0;
      if (char * stop = strpbrk(line+1, " \r\n"))
        headerlen = stop - (line+1);
      else
        headerlen = strlen(line+1);
      
      headerchars += headerlen;
      
      if (headerlen > longestheader)
        longestheader = headerlen;


      /* store the line number */
      
      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += MEMCHUNK;
          datap = (char *) xrealloc(datap, dataalloc);
        }
      memcpy(datap + datalen, & lineno, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


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
      if (!fgets(line, LINEALLOC, fp))
        line[0] = 0;
      lineno++;


      /* read and store sequence */

      unsigned long seqbegin = datalen;

      while (line[0] && (line[0] != '>'))
        {
          unsigned char c;
          char * p = line;
          while((c = *p++))
	    {
	    char m;
            if ((m = map_nt[(unsigned int)c]) >= 0)
              {
                while (datalen >= dataalloc)
                  {
                    dataalloc += MEMCHUNK;
                    datap = (char *) xrealloc(datap, dataalloc);
                  }
                
                *(datap+datalen) = m;
                datalen++;
              }
            else if ((c != 10) && (c != 13))
              {
                if ((c >= 32) && (c <= 126))
                  fprintf(stderr,
                          "\nIllegal character '%c' in sequence on line %u\n",
                          c,
                          lineno);
                else
                  fprintf(stderr,
                          "\nIllegal character (ascii no %d) in sequence on line %u\n",
                          c,
                          lineno);
                exit(1);
              }
	    }
          line[0] = 0;
          if (!fgets(line, LINEALLOC, fp))
            line[0] = 0;
          lineno++;
        }
      
      while (datalen >= dataalloc)
        {
          dataalloc += MEMCHUNK;
          datap = (char *) xrealloc(datap, dataalloc);
        }
      
      long length = datalen - seqbegin;

      if (length == 0)
        {
          fprintf(stderr, "\nError: Empty sequence found on line %u.\n\n", lineno-1);
          exit(1);
        }

      nucleotides += length;

      if (length > longest)
        longest = length;

      *(datap+datalen) = 0;
      datalen++;

      sequences++;
      
      if (is_regular)
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

  /* set up hash to check for unique sequences */

  unsigned long seqhashsize = 2 * sequences;

  seqinfo_t * * seqhashtable = 0;

  if (opt_differences > 0)
    {
      seqhashtable =
        (seqinfo_t **) xmalloc(seqhashsize * sizeof(seqinfo_t *));
      memset(seqhashtable, 0, seqhashsize * sizeof(seqinfo_t *));
    }

  /* create indices */

  seqindex = (seqinfo_t *) xmalloc(sequences * sizeof(seqinfo_t));
  seqinfo_t * seqindex_p = seqindex;

  regex_t db_regexp;
  regmatch_t pmatch[4];

  if (opt_usearch_abundance)
    {
      if (regcomp(&db_regexp, "(^|;)size=([0-9]+)(;|$)", REG_EXTENDED))
        fatal("Regular expression compilation failed");
    }
  else
    {
      if (regcomp(&db_regexp, "(_)([0-9]+)$", REG_EXTENDED))
        fatal("Regular expression compilation failed");
    }

  seqinfo_t * lastseq = 0;

  int presorted = 1;
  int missingabundance = 0;
  unsigned int missingabundance_lineno = 0;
  char * missingabundance_header = 0;

  char * p = datap;
  progress_init("Indexing database:", sequences);
  for(unsigned long i=0; i<sequences; i++)
    {
      /* get line number */
      unsigned int lineno = *((unsigned int*)p);
      p += sizeof(unsigned int);

      /* get header */
      seqindex_p->header = p;
      seqindex_p->headerlen = strlen(seqindex_p->header);
      p += seqindex_p->headerlen + 1;

      /* and sequence */
      seqindex_p->seq = p;
      seqindex_p->seqlen = strlen(p);
      p += seqindex_p->seqlen + 1;

      /* get amplicon abundance */
      if (!regexec(&db_regexp, seqindex_p->header, 4, pmatch, 0))
        {
          seqindex_p->abundance = atol(seqindex_p->header + pmatch[2].rm_so);
          seqindex_p->abundance_start = pmatch[0].rm_so;
          seqindex_p->abundance_end = pmatch[0].rm_eo;

          if (seqindex_p->abundance == 0)
            {
              fprintf(stderr,
                      "\nError: Illegal abundance value on line %u:\n%s\n"
                      "Abundance values should be positive integers.\n\n",
                      lineno,
                      seqindex_p->header);
              exit(1);
            }
        }
      else
        {
          seqindex_p->abundance_start = seqindex_p->headerlen;
          seqindex_p->abundance_end = seqindex_p->headerlen;
          seqindex_p->abundance = 0;
        }
      
      if (seqindex_p->abundance < 1)
        {
          if (opt_append_abundance)
            {
              seqindex_p->abundance = opt_append_abundance;
            }
          else
            {
              missingabundance++;
              if (missingabundance == 1)
                {
                  missingabundance_lineno = lineno;
                  missingabundance_header = seqindex_p->header;
                }
            }
        }

      if (seqindex_p->abundance_start == 0)
          fatal("Empty sequence identifier");


      /* check if the sequences are presorted by abundance and header */

      if (presorted && lastseq)
        {
          if (lastseq->abundance < seqindex_p->abundance)
            presorted = 0;
          else if (lastseq->abundance == seqindex_p->abundance)
            {
              if (strcmp(lastseq->header, seqindex_p->header) > 0)
                presorted = 0;
            }
        }

      lastseq = seqindex_p;


      /* check for duplicated identifiers using hash table */

      unsigned long hdrhash = HASH((unsigned char*)seqindex_p->header, seqindex_p->abundance_start);
      seqindex_p->hdrhash = hdrhash;
      unsigned long hdrhashindex = hdrhash % hdrhashsize;

      seqinfo_t * hdrfound = 0;
    
      while ((hdrfound = hdrhashtable[hdrhashindex]))
        {
          if ((hdrfound->hdrhash == hdrhash) &&
              (hdrfound->abundance_start == seqindex_p->abundance_start) &&
              (strncmp(hdrfound->header, seqindex_p->header, hdrfound->abundance_start) == 0))
            break;
          hdrhashindex = (hdrhashindex + 1) % hdrhashsize;
        }

      if (hdrfound)
        {
          duplicatedidentifiers++;
          fprintf(stderr, "\nError: Duplicated sequence identifier: %.*s\n\n",
                  seqindex_p->abundance_start,
                  seqindex_p->header);
          exit(1);
        }

      hdrhashtable[hdrhashindex] = seqindex_p;
    

      if (opt_differences > 0)
        {
          /* check for duplicated sequences using hash table */
          unsigned long seqhash = HASH((unsigned char*)seqindex_p->seq,
                                       seqindex_p->seqlen);
          seqindex_p->seqhash = seqhash;
          unsigned long seqhashindex = seqhash % seqhashsize;
          seqinfo_t * seqfound = 0;

          while ((seqfound = seqhashtable[seqhashindex]))
            {
              if ((seqfound->seqhash == seqhash) &&
                  (seqfound->seqlen == seqindex_p->seqlen) &&
                  (memcmp(seqfound->seq, seqindex_p->seq, seqfound->seqlen) == 0))
                break;
              seqhashindex = (seqhashindex + 1) % seqhashsize;
            }

          if (seqfound)
            duplicates_found++;
          else
            seqhashtable[seqhashindex] = seqindex_p;
        }

      seqindex_p++;
      progress_update(i);
    }
  progress_done();

  if (missingabundance)
    {
      fprintf(stderr,
              "\nError: Abundance annotations not found for %d sequences, starting on line %u.\n"
              ">%s\n"
              "Fasta headers must end with abundance annotations (_INT or ;size=INT).\n"
              "The -z option must be used if the abundance annotation is in the latter format.\n"
              "Abundance annotations can be produced by dereplicating the sequences.\n"
              "The header is defined as the string comprised between the \">\" symbol\n"
              "and the first space or the end of the line, whichever comes first.\n\n",
              missingabundance,
              missingabundance_lineno,
              missingabundance_header);
      exit(1);
    }

  if (duplicates_found)
    {
      fprintf(logfile,
              "WARNING: %lu duplicated sequences detected.\n"
              "Please consider dereplicating your data for optimal results.\n",
              duplicates_found);
    }

  if (!presorted)
    {
      progress_init("Abundance sorting:", 1);
      qsort(seqindex, sequences, sizeof(seqinfo_t), db_compare_abundance);
      progress_done();
    }

  regfree(&db_regexp);

  free(hdrhashtable);

  if (seqhashtable)
    {
      free(seqhashtable);
      seqhashtable = 0;
    }
}

void db_qgrams_init()
{
  qgrams = (qgramvector_t *) xmalloc(sequences * sizeof(qgramvector_t));

  seqinfo_t * seqindex_p = seqindex;
  progress_init("Find qgram vects: ", sequences);
  for(unsigned int i=0; i<sequences; i++)
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

void db_fprintseq(FILE * fp, int a, int width)
{
  char * seq = db_getsequence(a);
  int len = db_getsequencelen(a);
  char buffer[1025];
  char * buf;

  if (len < 1025)
    buf = buffer;
  else
    buf = (char*) xmalloc(len+1);

  for(int i=0; i<len; i++)
    buf[i] = sym_nt[(int)(seq[i])];
  buf[len] = 0;

  if (width < 1)
    fprintf(fp, "%.*s\n", (int)(len), buf);
  else
    {
      long rest = len;
      for(int i=0; i<len; i += width)
        {
          fprintf(fp, "%.*s\n", (int)(MIN(rest,width)), buf+i);
          rest -= width;
        }
    }

  if (len >= 1025)
    free(buf);
}
