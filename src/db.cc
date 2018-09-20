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
unsigned int longest = 0;
int longestheader = 0;

seqinfo_t * seqindex = 0;
char * datap = 0;
qgramvector_t * qgrams = 0;

static int missingabundance = 0;
static unsigned long missingabundance_lineno = 0;
static char * missingabundance_header = 0;

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

bool find_swarm_abundance(const char * header,
                          int * start,
                          int * end,
                          long * number)
{
  /*
    Identify the first occurence of the pattern (_)([0-9]+)$
    in the header string.
  */

  const char * digit_chars = "0123456789";

  if (!header)
    return false;

  const char * us = strrchr(header, '_');

  if (!us)
    return false;

  size_t digits = strspn(us + 1, digit_chars);

  if (us[1 + digits] == 0)
    {
      * start = us - header;
      * end = *start + digits;
      * number = atol(us + 1);
      return true;
    }
  else
    {
      * start = 0;
      * end = 0;
      * number = 0;
      return false;
    }
}

bool find_usearch_abundance(const char * header,
                            int * start,
                            int * end,
                            long * number)
{
  /*
    Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
    in the header string.
  */

  const char * attribute = "size=";
  const char * digit_chars = "0123456789";

  if ((! header) || (! attribute))
    return false;

  int hlen = strlen(header);
  int alen = strlen(attribute);

  int i = 0;

  while (i < hlen - alen)
    {
      char * r = (char *) strstr(header + i, attribute);

      /* no match */
      if (r == NULL)
        break;

      i = r - header;

      /* check for ';' in front */
      if ((i > 0) && (header[i-1] != ';'))
        {
          i += alen + 1;
          continue;
        }

      int digits = (int) strspn(header + i + alen, digit_chars);

      /* check for at least one digit */
      if (digits == 0)
        {
          i += alen + 1;
          continue;
        }

      /* check for ';' after */
      if ((i + alen + digits < hlen) && (header[i + alen + digits] != ';'))
        {
          i += alen + digits + 2;
          continue;
        }

      /* ok */
      * start = MAX(0, i - 1);
      * end   = MIN(i + alen + digits + 1, hlen);
      * number = atol(header + i + alen);
      return true;
    }
  * start = 0;
  * end = 0;
  * number = 0;
  return false;
}

void find_abundance(struct seqinfo_s * sp, unsigned long lineno)
{
  char * header = sp->header;

  /* read size/abundance annotation */
  long abundance = 0;
  int start = 0;
  int end = 0;
  long number = 0;

  if (opt_usearch_abundance)
    {
      /* (^|;)size=([0-9]+)(;|$) */

      if (find_usearch_abundance(header, & start, & end, & number))
        {
          if (number > 0)
            abundance = number;
          else
            {
              fprintf(stderr,
                      "\nError: Illegal abundance value on line %lu:\n%s\n"
                      "Abundance values should be positive integers.\n\n",
                      lineno,
                      header);
              exit(1);
            }
        }
    }
  else
    {
      /* (_)([0-9]+)$ */

      if (find_swarm_abundance(header, & start, & end, & number))
        {
          if (number > 0)
            abundance = number;
          else
            {
              fprintf(stderr,
                      "\nError: Illegal abundance value on line %lu:\n%s\n"
                      "Abundance values should be positive integers.\n\n",
                      lineno,
                      header);
              exit(1);
            }
        }
    }

  if (abundance == 0)
    {
      start = strlen(header);
      end = start;

      if (opt_append_abundance)
        abundance = opt_append_abundance;
      else
        {
          missingabundance++;
          if (missingabundance == 1)
            {
              missingabundance_lineno = lineno;
              missingabundance_header = header;
            }
        }
    }

  sp->abundance = abundance;
  sp->abundance_start = start;
  sp->abundance_end = end;
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

  progress_init("Reading sequences:", filesize);

  while(line[0])
    {
      /* read header */
      /* the header ends at a space, cr, lf or null character */

      if (line[0] != '>')
        fatal("Illegal header line in fasta file.");

      long headerlen = strcspn(line + 1, " \r\n");

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


      /* store a dummy sequence length */

      unsigned int length = 0;

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += MEMCHUNK;
          datap = (char *) xrealloc(datap, dataalloc);
        }
      unsigned long datalen_seqlen = datalen;
      memcpy(datap + datalen, & length, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* read and store sequence */

      unsigned long nt_buffer = 0;
      unsigned int nt_bufferlen = 0;
      const unsigned int nt_buffersize = 4 * sizeof(nt_buffer);

      while (line[0] && (line[0] != '>'))
        {
          char c;
          char * p = line;
          while((c = *p++))
	    {
              char m;
              if ((m = map_nt[(int)c]) >= 0)
                {
                  nt_buffer |= (((unsigned long)m)-1) << (2 * nt_bufferlen);
                  length++;
                  nt_bufferlen++;

                  if (nt_bufferlen == nt_buffersize)
                    {
                      while (datalen + sizeof(nt_buffer) > dataalloc)
                        {
                          dataalloc += MEMCHUNK;
                          datap = (char *) xrealloc(datap, dataalloc);
                        }

                      memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
                      datalen += sizeof(nt_buffer);

                      nt_bufferlen = 0;
                      nt_buffer = 0;
                    }
                }
              else if ((c != 10) && (c != 13))
                {
                  if ((c >= 32) && (c <= 126))
                    fprintf(stderr,
                            "Illegal character '%c' in sequence on line %u",
                            c,
                            lineno);
                  else
                    fprintf(stderr,
                            "Illegal character (ascii no %d) in sequence on line %u",
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

      /* fill in real length */

      memcpy(datap + datalen_seqlen, & length, sizeof(unsigned int));

      if (length == 0)
        {
          fprintf(stderr, "\nError: Empty sequence found on line %u.\n\n", lineno-1);
          exit(1);
        }

      nucleotides += length;

      if (length > longest)
        longest = length;


      /* save remaining padded unsigned long, if any */

      if (nt_bufferlen > 0)
        {
          while (datalen + sizeof(nt_buffer) > dataalloc)
            {
              dataalloc += MEMCHUNK;
              datap = (char *) xrealloc(datap, dataalloc);
            }

          memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
          datalen += sizeof(nt_buffer);

          nt_buffer = 0;
          nt_bufferlen = 0;
        }

      sequences++;

      if (is_regular)
        progress_update(ftell(fp));
    }
  progress_done();

  fclose(fp);

  /* init zobrist hashing */

  zobrist_init(longest + 2);  // add 2 for two insertions

  /* set up hash to check for unique headers */

  unsigned long hdrhashsize = 2 * sequences;

  seqinfo_t * * hdrhashtable =
    (seqinfo_t **) xmalloc(hdrhashsize * sizeof(seqinfo_t *));
  memset(hdrhashtable, 0, hdrhashsize * sizeof(seqinfo_t *));

  unsigned long duplicatedidentifiers = 0;

  /* set up hash to check for unique sequences */

  unsigned long seqhashsize = 2 * sequences;

  seqinfo_t * * seqhashtable = 0;

  if (opt_differences > 1)
    {
      seqhashtable =
        (seqinfo_t **) xmalloc(seqhashsize * sizeof(seqinfo_t *));
      memset(seqhashtable, 0, seqhashsize * sizeof(seqinfo_t *));
    }

  /* create indices */

  seqindex = (seqinfo_t *) xmalloc(sequences * sizeof(seqinfo_t));
  seqinfo_t * seqindex_p = seqindex;

  seqinfo_t * lastseq = 0;

  int presorted = 1;

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
      unsigned int seqlen = *((unsigned int*)p);
      seqindex_p->seqlen = seqlen;
      p += sizeof(unsigned int);
      seqindex_p->seq = p;
      p += nt_bytelength(seqlen);

      /* get amplicon abundance */
      find_abundance(seqindex_p, lineno);

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

      unsigned long hdrhash = HASH((unsigned char*)seqindex_p->header,
                                   seqindex_p->abundance_start);
      seqindex_p->hdrhash = hdrhash;
      unsigned long hdrhashindex = hdrhash % hdrhashsize;

      seqinfo_t * hdrfound = 0;

      while ((hdrfound = hdrhashtable[hdrhashindex]))
        {
          if ((hdrfound->hdrhash == hdrhash) &&
              (hdrfound->abundance_start == seqindex_p->abundance_start) &&
              (strncmp(hdrfound->header,
                       seqindex_p->header,
                       hdrfound->abundance_start) == 0))
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

      /* hash sequence */
      seqindex_p->seqhash = zobrist_hash((unsigned char*)seqindex_p->seq,
                                         seqindex_p->seqlen);

      if (opt_differences > 1)
        {
          /* Check for duplicated sequences using hash table, */
          /* but only for d>1. Handled internally for d=1.    */

          unsigned long seqhashindex = seqindex_p->seqhash % seqhashsize;
          seqinfo_t * seqfound = 0;

          while ((seqfound = seqhashtable[seqhashindex]))
            {
              if ((seqfound->seqhash == seqindex_p->seqhash) &&
                  (seqfound->seqlen == seqindex_p->seqlen) &&
                  (memcmp(seqfound->seq,
                          seqindex_p->seq,
                          nt_bytelength(seqindex_p->seqlen)) == 0))
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
              "\nError: Abundance annotations not found for %d sequences, starting on line %lu.\n"
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

unsigned long db_gethash(unsigned long seqno)
{
  return seqindex[seqno].seqhash;
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
    putchar(sym_nt[1+nt_extract(seq, i)]);
}

void db_free()
{
  zobrist_exit();

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
    buf[i] = sym_nt[1+nt_extract(seq, i)];
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
