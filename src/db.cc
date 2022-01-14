/*
    SWARM

    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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

constexpr unsigned int memchunk {1 << 20};  // 1 megabyte
constexpr unsigned int linealloc {2048};
constexpr long unsigned int n_chars {INT8_MAX + 1};  // 128 ascii chars
constexpr auto max_sequence_length {67108861};  // (2^26 - 3)
// for longer sequences, 'zobrist_tab_byte_base' is bigger than 8 x
// 2^32 (512 x max_sequence_length) and cannot be addressed with
// uint32 pointers, which leads to a segmentation fault

auto make_nt_map () -> std::array<signed char, n_chars> {
    // set the 128 ascii chars to '-1' except Aa, Cc, Gg, Tt and Uu
  std::array<signed char, n_chars> ascii_map {{0}};
    ascii_map.fill(-1);
    ascii_map['A'] = 1;
    ascii_map['a'] = 1;
    ascii_map['C'] = 2;
    ascii_map['c'] = 2;
    ascii_map['G'] = 3;
    ascii_map['g'] = 3;
    ascii_map['T'] = 4;
    ascii_map['t'] = 4;
    ascii_map['U'] = 4;
    ascii_map['u'] = 4;
    return ascii_map;
    }

const auto map_nt = make_nt_map();

const std::array<char, 32> sym_nt =
  {'-', 'A', 'C', 'G', 'T', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
   ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};


static unsigned int sequences {0};
static uint64_t nucleotides {0};
static uint64_t headerchars {0};
static unsigned int longest {0};
static unsigned int longestheader {0};
static char * datap {nullptr};
static int missingabundance {0};
static uint64_t missingabundance_lineno {0};
static char * missingabundance_header {nullptr};

struct seqinfo_s
{
  char * header;
  char * seq;
  uint64_t abundance;
  uint64_t hdrhash;
  uint64_t seqhash;
  int headerlen;
  unsigned int seqlen;
  unsigned int clusterid;
  int abundance_start;
  int abundance_end;
  int dummy; /* alignment padding only */
};

using seqinfo_t = struct seqinfo_s;
extern seqinfo_t * seqindex;

seqinfo_t * seqindex {nullptr};
qgramvector_t * qgrams {nullptr};


auto db_getnucleotidecount() -> uint64_t
{
  return nucleotides;
}


auto db_getsequencecount() -> unsigned int
{
  return sequences;
}


auto db_getlongestsequence() -> unsigned int
{
  return longest;
}


auto fprint_id(std::FILE * stream, uint64_t x, bool opt_usearch_abundance,
               int64_t opt_append_abundance) -> void
{
  const seqinfo_t * sp = seqindex + x;
  const char * h = sp->header;
  const int hdrlen = sp->headerlen;

  if ((opt_append_abundance != 0) && (sp->abundance_start == sp->abundance_end)) {
    if (opt_usearch_abundance) {
      fprintf(stream, "%.*s;size=%" PRIu64 ";", hdrlen, h, sp->abundance);
    }
    else {
      fprintf(stream, "%.*s_%" PRIu64, hdrlen, h, sp->abundance);
    }
  }
  else {
    fprintf(stream, "%.*s", hdrlen, h);
  }
}


auto fprint_id_noabundance(std::FILE * stream, uint64_t seqno, bool opt_usearch_abundance) -> void
{
  const seqinfo_t * sp = seqindex + seqno;
  const char * h = sp->header;
  const int hdrlen = sp->headerlen;

  if (sp->abundance_start < sp->abundance_end)
    {
      /* print start of header */
      fprintf(stream, "%.*s", sp->abundance_start, h);

      if (opt_usearch_abundance)
        {
          /* print semicolon if the abundance is not at either end */
          if ((sp->abundance_start > 0) && (sp->abundance_end < hdrlen)) {
            fprintf(stream, ";");
          }

          /* print remaining part */
          fprintf(stream, "%.*s", hdrlen - sp->abundance_end, h + sp->abundance_end);
        }
    }
  else {
    fprintf(stream, "%.*s", hdrlen, h);
  }
}


auto fprint_id_with_new_abundance(std::FILE * stream,
                                  uint64_t seqno,
                                  uint64_t abundance,
                                  bool opt_usearch_abundance) -> void
{
  const seqinfo_t * sp = seqindex + seqno;

  if (opt_usearch_abundance) {
    fprintf(stream,
            "%.*s%ssize=%" PRIu64 ";%.*s",
            sp->abundance_start,
            sp->header,
            sp->abundance_start > 0 ? ";" : "",
            abundance,
            sp->headerlen - sp->abundance_end,
            sp->header + sp->abundance_end);
  }
  else {
    fprintf(stream,
            "%.*s_%" PRIu64,
            sp->abundance_start,
            sp->header,
            abundance);
  }
}


auto db_compare_abundance(const void * a, const void * b) -> int
{
  const auto * x = reinterpret_cast<const seqinfo_t *>(a);
  const auto * y = reinterpret_cast<const seqinfo_t *>(b);
  int status {0};

  if (x->abundance > y->abundance) {
    status = -1;
  }
  else if (x->abundance < y->abundance) {
    status = +1;
  }
  else {
    status = strcmp(x->header, y->header);
  }
  return status;
}


auto find_swarm_abundance(const char * header,
                          int * start,
                          int * end,
                          int64_t * number) -> bool
{
  /*
    Identify the first occurence of the pattern (_)([0-9]+)$
    in the header string.
  */

  * start = 0;
  * end = 0;
  * number = 0;

  constexpr unsigned int max_digits {20};  // 20 digits at most (abundance > 10^20)
  const std::string digit_chars = "0123456789";

  if (header == nullptr) {
    return false;
  }

  const char * us = strrchr(header, '_');

  if (us == nullptr) {
    return false;
  }

  size_t digits = strspn(us + 1, digit_chars.c_str());

  if (digits > max_digits) {
    return false;
  }

  if (us[digits + 1] != 0) {
    return false;
  }

  int64_t s = us - header;
  int64_t e = s + 1 + static_cast<int64_t>(digits);

  * start = static_cast<int>(s);
  * end = static_cast<int>(e);
  * number = atol(us + 1);

  return true;
}


auto find_usearch_abundance(const char * header,
                            int * start,
                            int * end,
                            int64_t * number) -> bool
{
  /*
    Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
    in the header string.
  */

  assert(header != nullptr); // header cannot be a nullptr at this stage

  const std::string attribute {"size="};
  const std::string digit_chars {"0123456789"};
  const uint64_t hlen = strlen(header);
  const uint64_t alen = attribute.length();
  uint64_t i = 0;

  while (i + alen < hlen)
    {
      const char * r = strstr(header + i, attribute.c_str());

      /* no match */
      if (r == nullptr) {
        break;
      }

      i = static_cast<uint64_t>(r - header);

      /* check for ';' in front */
      if ((i > 0) && (header[i-1] != ';'))
        {
          i += alen + 1;
          continue;
        }

      uint64_t digits = strspn(header + i + alen, digit_chars.c_str());

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
      if (i > 0) {
        * start = static_cast<int>(i - 1);
      }
      else {
        * start = 0;
      }
      * end   = static_cast<int>(std::min(i + alen + digits + 1, hlen));
      * number = atol(header + i + alen);
      return true;
    }

  return false;
}


auto msg_illegal_abundance_value(uint64_t lineno, char * header) -> void {
  fprintf(stderr,
          "\nError: Illegal abundance value on line %" PRIu64 ":\n%s\n"
          "Abundance values should be positive integers.\n\n",
          lineno,
          header);
}


void find_abundance(struct seqinfo_s * sp, uint64_t lineno, bool opt_usearch_abundance, int64_t opt_append_abundance)
{
  char * header = sp->header;

  /* read size/abundance annotation */
  int64_t abundance = 0;
  int start = 0;
  int end = 0;
  int64_t number = 0;

  if (opt_usearch_abundance)
    {
      /* (^|;)size=([0-9]+)(;|$) */

      if (find_usearch_abundance(header, & start, & end, & number))
        {
          if (number <= 0) {
            msg_illegal_abundance_value(lineno, header);
            exit(1);
          }
          abundance = number;
        }
    }
  else
    {
      /* (_)([0-9]+)$ */

      if (find_swarm_abundance(header, & start, & end, & number))
        {
          if (number <= 0) {
            msg_illegal_abundance_value(lineno, header);
            exit(1);
          }
          abundance = number;
        }
    }

  if (abundance == 0)
    {
      start = sp->headerlen;
      end = start;

      if (opt_append_abundance != 0) {
        abundance = opt_append_abundance;
      }
      else
        {
          missingabundance++;
          // record the position of the first missing abundance entry
          if (missingabundance == 1)
            {
              missingabundance_lineno = lineno;
              missingabundance_header = header;
            }
        }
    }

  sp->abundance = static_cast<uint64_t>(abundance);
  sp->abundance_start = start;
  sp->abundance_end = end;
}


void db_read(const char * filename, struct Parameters const & p)
{
  /* allocate space */

  uint64_t dataalloc {memchunk};
  datap = static_cast<char *>(xmalloc(dataalloc));
  uint64_t datalen {0};
  uint64_t duplicates_found {0};

  longest = 0;
  longestheader = 0;
  sequences = 0;
  nucleotides = 0;
  headerchars = 0;

  std::FILE * fp {nullptr};
  if (filename != nullptr)
    {
      fp = fopen_input(filename);
      if (fp == nullptr)
        {
          fprintf(stderr, "\nError: Unable to open input data file (%s).\n", filename);
          exit(1);
        }
    }
  else {
    fp = stdin;
  }

  /* get file size */

  struct stat fs;

  if (fstat(fileno(fp), & fs) != 0)
    {
      fprintf(stderr, "\nUnable to fstat on input file (%s)\n", filename);
      exit(1);
    }
  bool is_regular = S_ISREG(fs.st_mode);
  uint64_t filesize = is_regular ? fs.st_size : 0;
  uint64_t filepos = 0;

  if (! is_regular) {
    fprintf(logfile, "Waiting for data... (Hit Ctrl-C and run swarm -h if you meant to read data from a file.)\n");
  }

  size_t linecap = linealloc;
  char * line = static_cast<char*>(xmalloc(linecap));
  ssize_t linelen = xgetline(& line, & linecap, fp);
  if (linelen < 0)
    {
      line[0] = 0;
      linelen = 0;
    }
  filepos += linelen;

  unsigned int lineno {1};

  progress_init("Reading sequences:", filesize);

  while(line[0] != 0)
    {
      /* read header */
      /* the header ends at a space, cr, lf or null character */

      if (line[0] != '>') {
        fatal("Illegal header line in fasta file.");
      }

      auto headerlen = static_cast<unsigned int>
        (strcspn(line + 1, " \r\n"));

      headerchars += headerlen;

      if (headerlen > longestheader) {
        longestheader = headerlen;
      }


      /* store the line number */

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      memcpy(datap + datalen, & lineno, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* store the header */

      while (datalen + headerlen + 1 > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      memcpy(datap + datalen, line + 1, headerlen);
      *(datap + datalen + headerlen) = 0;
      datalen += headerlen + 1;


      /* get next line */

      linelen = xgetline(& line, & linecap, fp);
      if (linelen < 0)
        {
          line[0] = 0;
          linelen = 0;
        }
      filepos += linelen;

      lineno++;


      /* store a dummy sequence length */

      unsigned int length {0};

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      uint64_t datalen_seqlen = datalen;
      memcpy(datap + datalen, & length, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* read and store sequence */

      uint64_t nt_buffer {0};
      unsigned int nt_bufferlen {0};
      const unsigned int nt_buffersize {4 * sizeof(nt_buffer)};
      constexpr int new_line {10};
      constexpr int carriage_return {13};
      constexpr int start_chars_range {32};  // visible ascii chars: 32-126
      constexpr int end_chars_range {126};

      while ((line[0] != 0) && (line[0] != '>'))
        {
          unsigned char c {0};
          char * pl = line;
          while((c = static_cast<unsigned char>(*pl++)) != 0U)
            {
              signed char m {0};
              if ((m = map_nt[static_cast<unsigned int>(c)]) >= 0)
                {
                  nt_buffer |= ((static_cast<uint64_t>(m))-1) << (2 * nt_bufferlen);
                  length++;
                  nt_bufferlen++;

                  if (nt_bufferlen == nt_buffersize)
                    {
                      while (datalen + sizeof(nt_buffer) > dataalloc)
                        {
                          dataalloc += memchunk;
                          datap = static_cast<char *>(xrealloc(datap, dataalloc));
                        }

                      memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
                      datalen += sizeof(nt_buffer);

                      nt_bufferlen = 0;
                      nt_buffer = 0;
                    }
                }
              else if ((c != new_line) && (c != carriage_return))
                {
                  if ((c >= start_chars_range) && (c <= end_chars_range)) {
                    fprintf(stderr,
                            "\nError: Illegal character '%c' in sequence on line %u\n",
                            c,
                            lineno);
                  }
                  else {
                    fprintf(stderr,
                            "\nError: Illegal character (ascii no %d) in sequence on line %u\n",
                            c,
                            lineno);
                  }
                  exit(1);
                }
            }

          /* check length of longest sequence */
          if (length > max_sequence_length) {
            fatal("Sequences longer than 67 108 861 symbols are not supported.");
          }

          linelen = xgetline(& line, & linecap, fp);
          if (linelen < 0)
            {
              line[0] = 0;
              linelen = 0;
            }
          filepos += linelen;

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

      if (length > longest) {
        longest = length;
      }


      /* save remaining padded 64-bit value with nt's, if any */

      if (nt_bufferlen > 0)
        {
          while (datalen + sizeof(nt_buffer) > dataalloc)
            {
              dataalloc += memchunk;
              datap = static_cast<char *>(xrealloc(datap, dataalloc));
            }

          memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
          datalen += sizeof(nt_buffer);

          nt_buffer = 0;
          nt_bufferlen = 0;  // that value is never read again, all tests pass without it
        }

      sequences++;

      if (is_regular) {
        progress_update(filepos);
      }
    }
  progress_done();

  fclose(fp);

  /* init zobrist hashing */

  // add 2 for two insertions
  unsigned int zobrist_len = std::max(4 * longestheader, longest + 2);
  zobrist_init(zobrist_len);

  /* set up hash to check for unique headers */

  uint64_t hdrhashsize {2 * sequences};

  auto * * hdrhashtable =
    static_cast<seqinfo_t **>(xmalloc(hdrhashsize * sizeof(seqinfo_t *)));
  memset(hdrhashtable, 0, hdrhashsize * sizeof(seqinfo_t *));

  uint64_t duplicatedidentifiers {0};

  /* set up hash to check for unique sequences */

  uint64_t seqhashsize {2 * sequences};

  seqinfo_t * * seqhashtable {nullptr};

  if (p.opt_differences > 1)
    {
      seqhashtable =
        static_cast<seqinfo_t **>(xmalloc(seqhashsize * sizeof(seqinfo_t *)));
      memset(seqhashtable, 0, seqhashsize * sizeof(seqinfo_t *));
    }

  /* create indices */

  seqindex = static_cast<seqinfo_t *>(xmalloc(sequences * sizeof(seqinfo_t)));
  seqinfo_t * seqindex_p {seqindex};

  seqinfo_t * lastseq {nullptr};

  bool presorted {true};

  char * pl {datap};
  progress_init("Indexing database:", sequences);
  for(auto i = 0ULL; i < sequences; i++)
    {
      /* get line number */
      unsigned int line_number = *(reinterpret_cast<unsigned int*>(pl));
      pl += sizeof(unsigned int);

      /* get header */
      seqindex_p->header = pl;
      seqindex_p->headerlen = static_cast<int>(strlen(seqindex_p->header));
      pl += seqindex_p->headerlen + 1;

      /* and sequence */
      unsigned int seqlen = *(reinterpret_cast<unsigned int*>(pl));
      seqindex_p->seqlen = seqlen;
      pl += sizeof(unsigned int);
      seqindex_p->seq = pl;
      pl += nt_bytelength(seqlen);

      /* get amplicon abundance */
      find_abundance(seqindex_p, line_number, p.opt_usearch_abundance, p.opt_append_abundance);

      if ((seqindex_p->abundance_start == 0) &&
          (seqindex_p->abundance_end == seqindex_p->headerlen)) {
        fatal("Empty sequence identifier");
      }

      /* check if the sequences are presorted by abundance and header */

      if (presorted && (lastseq != nullptr))
        {
          if (lastseq->abundance < seqindex_p->abundance) {
            presorted = false;
          }
          else if (lastseq->abundance == seqindex_p->abundance)
            {
              if (strcmp(lastseq->header, seqindex_p->header) > 0) {
                presorted = false;
              }
            }
        }

      lastseq = seqindex_p;

      /* check for duplicated identifiers using hash table */

      /* find position and length of identifier in header */

      int id_start {0};
      int id_len {0};

      if (seqindex_p->abundance_start > 0)
        {
          /* id first, then abundance (e.g. >name;size=1 or >name_1) */
          id_start = 0;
          id_len = seqindex_p->abundance_start;
        }
      else
        {
          /* abundance first then id (e.g. >size=1;name) */
          id_start = seqindex_p->abundance_end;
          id_len = seqindex_p->headerlen - seqindex_p->abundance_end;
        }

      uint64_t hdrhash = zobrist_hash(reinterpret_cast<unsigned char*>
                                      (seqindex_p->header + id_start),
                                      4 * static_cast<unsigned int>(id_len));

      seqindex_p->hdrhash = hdrhash;
      uint64_t hdrhashindex = hdrhash % hdrhashsize;

      seqinfo_t * hdrfound {nullptr};

      while ((hdrfound = hdrhashtable[hdrhashindex]) != nullptr)
        {
          if (hdrfound->hdrhash == hdrhash)
            {
              int hit_id_start {0};
              int hit_id_len {0};

              if (hdrfound->abundance_start > 0)
                {
                  hit_id_start = 0;
                  hit_id_len = hdrfound->abundance_start;
                }
              else
                {
                  hit_id_start = hdrfound->abundance_end;
                  hit_id_len = hdrfound->headerlen - hdrfound->abundance_end;
                }

              if ((id_len == hit_id_len) &&
                  (strncmp(seqindex_p->header + id_start,
                           hdrfound->header + hit_id_start,
                           static_cast<uint64_t>(id_len)) == 0)) {
                break;
              }
            }

          hdrhashindex = (hdrhashindex + 1) % hdrhashsize;
        }

      if (hdrfound != nullptr)
        {
          duplicatedidentifiers++;
          fprintf(stderr, "\nError: Duplicated sequence identifier: %.*s\n\n",
                  id_len,
                  seqindex_p->header + id_start);
          exit(1);
        }

      hdrhashtable[hdrhashindex] = seqindex_p;

      /* hash sequence */
      seqindex_p->seqhash = zobrist_hash(reinterpret_cast<unsigned char*>
                                         (seqindex_p->seq),
                                         seqindex_p->seqlen);

      if (p.opt_differences > 1)
        {
          /* Check for duplicated sequences using hash table, */
          /* but only for d>1. Handled internally for d=1.    */

          uint64_t seqhashindex = seqindex_p->seqhash % seqhashsize;
          seqinfo_t * seqfound {nullptr};

          while ((seqfound = seqhashtable[seqhashindex]) != nullptr)
            {
              if ((seqfound->seqhash == seqindex_p->seqhash) &&
                  (seqfound->seqlen == seqindex_p->seqlen) &&
                  (memcmp(seqfound->seq,
                          seqindex_p->seq,
                          nt_bytelength(seqindex_p->seqlen)) == 0)) {
                break;
              }
              seqhashindex = (seqhashindex + 1) % seqhashsize;
            }

          if (seqfound != nullptr)
            {
              duplicates_found++;
              break;
            }
          seqhashtable[seqhashindex] = seqindex_p;
        }

      seqindex_p++;
      progress_update(i);
    }

  if (duplicates_found != 0U)
    {
      fprintf(logfile,
              "\n\n"
              "Error: some fasta entries have identical sequences.\n"
              "Swarm expects dereplicated fasta files.\n"
              "Such files can be produced with swarm or vsearch:\n"
              " swarm -d 0 -w derep.fasta -o /dev/null input.fasta\n"
              "or\n"
              " vsearch --derep_fulllength input.fasta --sizein --sizeout --output derep.fasta\n");
      exit(1);
    }

  progress_done();

  if (line != nullptr)
    {
      xfree(line);
      line = nullptr;
      linecap = 0;
    }

  if (missingabundance != 0)
    {
      fprintf(stderr,
              "\nError: Abundance annotations not found for %d sequences, starting on line %" PRIu64 ".\n"
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

  if (!presorted)
    {
      progress_init("Abundance sorting:", 1);
      qsort(seqindex, sequences, sizeof(seqinfo_t), db_compare_abundance);
      progress_done();
    }

  xfree(hdrhashtable);

  if (seqhashtable != nullptr)
    {
      xfree(seqhashtable);
      seqhashtable = nullptr;
    }

  // user report
  fprintf(logfile, "Database info:     %" PRIu64 " nt", db_getnucleotidecount());
  fprintf(logfile, " in %u sequences,", db_getsequencecount());
  fprintf(logfile, " longest %u nt\n", db_getlongestsequence());
}


void db_qgrams_init()
{
  qgrams = static_cast<qgramvector_t *>
    (xmalloc(sequences * sizeof(qgramvector_t)));

  seqinfo_t * seqindex_p {seqindex};
  progress_init("Find qgram vects: ", sequences);
  for(auto i = 0U; i < sequences; i++)
    {
      /* find qgrams */
      findqgrams(reinterpret_cast<unsigned char*>(seqindex_p->seq),
                 seqindex_p->seqlen,
                 qgrams[i]);
      seqindex_p++;
      progress_update(i);
    }
  progress_done();
}


void db_qgrams_done()
{
  xfree(qgrams);
}


auto db_gethash(uint64_t seqno) -> uint64_t
{
  return seqindex[seqno].seqhash;
}


auto db_getsequence(uint64_t seqno) -> char *
{
  return seqindex[seqno].seq;
}


void db_getsequenceandlength(uint64_t seqno,
                             char ** address,
                             unsigned int * length)
{
  *address = seqindex[seqno].seq;
  *length = seqindex[seqno].seqlen;
}


auto db_getsequencelen(uint64_t seqno) -> unsigned int
{
  return seqindex[seqno].seqlen;
}


auto db_getheader(uint64_t seqno) -> char *
{
  return seqindex[seqno].header;
}


auto db_getabundance(uint64_t seqno) -> uint64_t
{
  return seqindex[seqno].abundance;
}


void db_free()
{
  zobrist_exit();

  if (datap != nullptr) {
    xfree(datap);
  }
  if (seqindex != nullptr) {
    xfree(seqindex);
  }
}


auto db_fprintseq(std::FILE * fp, unsigned int a, unsigned int width) -> void
{
  constexpr static int default_length {1025};
  const unsigned int len = db_getsequencelen(a);
  char buffer[default_length];
  char * buf {nullptr};

  // buf = len < default_length ? buffer : static_cast<char*>(xmalloc(len+1));
  if (len < default_length) {
    buf = buffer;
  }
  else {
    buf = static_cast<char*>(xmalloc(len + 1));
  }

  for(auto i = 0U; i < len; i++) {
    buf[i] = sym_nt[1 + nt_extract(db_getsequence(a), i)];
  }
  buf[len] = 0;

  if (width < 1) {
    fprintf(fp, "%.*s\n", len, buf);
  }
  else { // unreachable: 'width' is always set to zero in algo, algod1 and derep
    auto rest = len;
    for(auto i = 0U; i < len; i += width) {
      fprintf(fp, "%.*s\n", std::min(rest, width), buf + i);
      rest -= width;
    }
  }

  if (len >= default_length) {
    xfree(buf);
  }
}
