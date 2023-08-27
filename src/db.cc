/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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
#include "qgram.h"
#include "util.h"
#include "utils/fatal.h"
#include "utils/nt_codec.h"
#include "zobrist.h"
#include <algorithm>  // std::min()
#include <array>
#include <cassert>  // assert()
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fclose()
#include <cstdlib>  // qsort()
#include <string>
#include <vector>


constexpr unsigned int memchunk {1U << 20U};  // 1 megabyte
constexpr unsigned int linealloc {2048};
constexpr long unsigned int n_chars {INT8_MAX + 1};  // 128 ascii chars
constexpr unsigned int max_sequence_length {67108861};  // (2^26 - 3)
// for longer sequences, 'zobrist_tab_byte_base' is bigger than 8 x
// 2^32 (512 x max_sequence_length) and cannot be addressed with
// uint32 pointers, which leads to a segmentation fault
constexpr unsigned int max_header_length {16777216 - 1};  // 2^24 minus 1

auto make_nt_map () -> std::array<uint64_t, n_chars> {
    // set the 128 ascii chars to zero except Aa, Cc, Gg, Tt and Uu
  std::array<uint64_t, n_chars> ascii_map {{0}};
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


auto fprint_id(std::FILE * stream, uint64_t seqno, bool opt_usearch_abundance,
               int64_t opt_append_abundance) -> void
{
  const seqinfo_t * seqinfo = seqindex + seqno;
  const char * hdrstr = seqinfo->header;
  const int hdrlen = seqinfo->headerlen;

  if ((opt_append_abundance != 0) && (seqinfo->abundance_start == seqinfo->abundance_end)) {
    if (opt_usearch_abundance) {
      std::fprintf(stream, "%.*s;size=%" PRIu64 ";", hdrlen, hdrstr, seqinfo->abundance);
    }
    else {
      std::fprintf(stream, "%.*s_%" PRIu64, hdrlen, hdrstr, seqinfo->abundance);
    }
  }
  else {
    std::fprintf(stream, "%.*s", hdrlen, hdrstr);
  }
}


auto fprint_id_noabundance(std::FILE * stream, uint64_t seqno, bool opt_usearch_abundance) -> void
{
  const seqinfo_t * seqinfo = seqindex + seqno;
  const char * hdrstr = seqinfo->header;
  const int hdrlen = seqinfo->headerlen;

  if (seqinfo->abundance_start < seqinfo->abundance_end)
    {
      /* print start of header */
      std::fprintf(stream, "%.*s", seqinfo->abundance_start, hdrstr);

      if (opt_usearch_abundance)
        {
          /* print semicolon if the abundance is not at either end */
          if ((seqinfo->abundance_start > 0) && (seqinfo->abundance_end < hdrlen)) {
            std::fprintf(stream, ";");
          }

          /* print remaining part */
          std::fprintf(stream, "%.*s", hdrlen - seqinfo->abundance_end, hdrstr + seqinfo->abundance_end);
        }
    }
  else {
    std::fprintf(stream, "%.*s", hdrlen, hdrstr);
  }
}


auto fprint_id_with_new_abundance(std::FILE * stream,
                                  uint64_t seqno,
                                  uint64_t abundance,
                                  bool opt_usearch_abundance) -> void
{
  const seqinfo_t * seqinfo = seqindex + seqno;

  if (opt_usearch_abundance) {
    std::fprintf(stream,
            "%.*s%ssize=%" PRIu64 ";%.*s",
            seqinfo->abundance_start,
            seqinfo->header,
            seqinfo->abundance_start > 0 ? ";" : "",
            abundance,
            seqinfo->headerlen - seqinfo->abundance_end,
            seqinfo->header + seqinfo->abundance_end);
  }
  else {
    std::fprintf(stream,
            "%.*s_%" PRIu64,
            seqinfo->abundance_start,
            seqinfo->header,
            abundance);
  }
}


auto db_compare_abundance(const void * ptr_a, const void * ptr_b) -> int
{
  const auto * seqinfo_x = reinterpret_cast<const seqinfo_t *>(ptr_a);
  const auto * seqinfo_y = reinterpret_cast<const seqinfo_t *>(ptr_b);

  if (seqinfo_x->abundance > seqinfo_y->abundance) {
    return -1;
  }

  if (seqinfo_x->abundance < seqinfo_y->abundance) {
    return +1;
  }

  return std::strcmp(seqinfo_x->header, seqinfo_y->header);
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

  static constexpr unsigned int max_digits {20};  // 20 digits at most (abundance > 10^20)
  static const std::string digit_chars = "0123456789";

  if (header == nullptr) {
    return false;  // refactoring: if header cannot be a nullptr, replace with assert
  }

  const char * abundance_string = std::strrchr(header, '_');

  if (abundance_string == nullptr) {
    return false;
  }

  const size_t n_digits = std::strspn(abundance_string + 1, digit_chars.c_str());

  if (n_digits > max_digits) {
    return false;
  }

  if (abundance_string[n_digits + 1] != 0) {
    return false;
  }

  const int64_t abundance_start = abundance_string - header;
  const int64_t abundance_end = abundance_start + 1 + static_cast<int64_t>(n_digits);

  * start = static_cast<int>(abundance_start);
  * end = static_cast<int>(abundance_end);
  * number = std::atol(abundance_string + 1);

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

  static const std::string attribute {"size="};
  static const std::string digit_chars {"0123456789"};
  const uint64_t hlen = strlen(header);
  const uint64_t alen = attribute.length();
  uint64_t position = 0;

  while (position + alen < hlen)
    {
      const char * result = std::strstr(header + position, attribute.c_str());

      /* no match */
      if (result == nullptr) {
        break;
      }

      position = static_cast<uint64_t>(result - header);

      /* check for ';' in front */
      if ((position > 0) && (header[position - 1] != ';'))
        {
          position += alen + 1;
          continue;
        }

      const uint64_t n_digits = std::strspn(header + position + alen, digit_chars.c_str());

      /* check for at least one digit */
      if (n_digits == 0)
        {
          position += alen + 1;
          continue;
        }

      /* check for ';' after */
      if ((position + alen + n_digits < hlen) && (header[position + alen + n_digits] != ';'))
        {
          position += alen + n_digits + 2;
          continue;
        }

      /* ok */
      if (position > 0) {
        * start = static_cast<int>(position - 1);
      }
      else {
        * start = 0;
      }
      * end   = static_cast<int>(std::min(position + alen + n_digits + 1, hlen));
      * number = std::atol(header + position + alen);

      return true;
    }

  return false;
}


void find_abundance(struct seqinfo_s * seqinfo, uint64_t lineno,
                    bool opt_usearch_abundance, int64_t opt_append_abundance)
{
  char * header = seqinfo->header;

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
            fatal(error_prefix, "Illegal abundance value on line ", lineno, ":\n",
                  header, "\nAbundance values should be positive integers.");
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
            fatal(error_prefix, "Illegal abundance value on line ", lineno, ":\n",
                  header, "\nAbundance values should be positive integers.");
          }
          abundance = number;
        }
    }

  if (abundance == 0)
    {
      start = seqinfo->headerlen;
      end = start;

      if (opt_append_abundance != 0) {
        abundance = opt_append_abundance;
      }
      else
        {
          ++missingabundance;
          // record the position of the first missing abundance entry
          if (missingabundance == 1)
            {
              missingabundance_lineno = lineno;
              missingabundance_header = header;
            }
        }
    }

  seqinfo->abundance = static_cast<uint64_t>(abundance);
  seqinfo->abundance_start = start;
  seqinfo->abundance_end = end;
}


auto db_read(const char * filename, struct Parameters const & parameters) -> void
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

  /* open input file or stream */

  assert(filename != nullptr);  // filename is set to '-' (stdin) by default

  std::FILE * input_fp { fopen_input(filename) };
  if (input_fp == nullptr)
    {
      fatal(error_prefix, "Unable to open input data file (", filename, ").\n");
    }

  /* get file size */

  struct stat fstat_buffer;  // refactor: complex uninitialized struct for fstat() 

  if (fstat(fileno(input_fp), &fstat_buffer) != 0)  // refactor: fstat and fileno linuxisms
    {
      fatal(error_prefix, "Unable to fstat on input file (", filename, ").\n");
    }
  const bool is_regular = S_ISREG(fstat_buffer.st_mode);  // refactor: S_ISREG linuxisms
  const uint64_t filesize = is_regular ? fstat_buffer.st_size : 0;
  uint64_t filepos = 0;

  if (! is_regular) {
    std::fprintf(logfile, "Waiting for data... (hit Ctrl-C and run 'swarm -h' if you meant to read data from a file)\n");
  }

  size_t linecap = linealloc;
  char * line {new char[linecap]};
  ssize_t linelen = xgetline(& line, & linecap, input_fp);
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
        fatal(error_prefix, "Illegal header line in fasta file.");
      }

      auto headerlen = static_cast<unsigned int>
        (std::strcspn(line + 1, " \r\n"));

      headerchars += headerlen;

      if (headerlen > longestheader) {
        longestheader = headerlen;
      }

      if (longestheader > max_header_length) {
        fatal(error_prefix, "Headers longer than 16,777,215 symbols are not supported.");
      }

      /* store the line number */

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      std::memcpy(datap + datalen, & lineno, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* store the header */

      while (datalen + headerlen + 1 > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      std::memcpy(datap + datalen, line + 1, headerlen);
      *(datap + datalen + headerlen) = 0;
      datalen += headerlen + 1;


      /* get next line */

      linelen = xgetline(& line, & linecap, input_fp);
      if (linelen < 0)
        {
          line[0] = 0;
          linelen = 0;
        }
      filepos += linelen;

      ++lineno;


      /* store a dummy sequence length */

      unsigned int length {0};

      while (datalen + sizeof(unsigned int) > dataalloc)
        {
          dataalloc += memchunk;
          datap = static_cast<char *>(xrealloc(datap, dataalloc));
        }
      const uint64_t datalen_seqlen = datalen;
      std::memcpy(datap + datalen, & length, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* read and store sequence */

      uint64_t nt_buffer {0};
      unsigned int nt_bufferlen {0};
      const unsigned int nt_buffersize {4 * sizeof(nt_buffer)};
      static constexpr int new_line {10};
      static constexpr int carriage_return {13};
      static constexpr int start_chars_range {32};  // visible ascii chars: 32-126
      static constexpr int end_chars_range {126};

      while ((line[0] != 0) && (line[0] != '>'))
        {
          unsigned char character {0};
          char * line_ptr = line;
          while((character = static_cast<signed char>(*line_ptr++)) != 0)
            {
              const auto mapped_char = map_nt[character];
              if (mapped_char != 0)
                {
                  nt_buffer |= (mapped_char - 1) << (2 * nt_bufferlen);
                  ++length;
                  ++nt_bufferlen;

                  if (nt_bufferlen == nt_buffersize)
                    {
                      while (datalen + sizeof(nt_buffer) > dataalloc)
                        {
                          dataalloc += memchunk;
                          datap = static_cast<char *>(xrealloc(datap, dataalloc));
                        }

                      std::memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
                      datalen += sizeof(nt_buffer);

                      nt_bufferlen = 0;
                      nt_buffer = 0;
                    }
                }
              else if ((character != new_line) && (character != carriage_return))
                {
                  if ((character >= start_chars_range) && (character <= end_chars_range)) {
                    fatal(error_prefix, "Illegal character '", character,
                          "' in sequence on line ", lineno, ".");
                  }
                  else {
                    fatal(error_prefix, "Illegal character (ascii no ", character,
                          ") in sequence on line ", lineno, ".");
                  }
                }
            }

          /* check length of longest sequence */
          if (length > max_sequence_length) {
            fatal(error_prefix, "Sequences longer than 67,108,861 symbols are not supported.");
          }

          linelen = xgetline(& line, & linecap, input_fp);
          if (linelen < 0)
            {
              line[0] = 0;
              linelen = 0;
            }
          filepos += linelen;

          ++lineno;
        }

      /* fill in real length */

      std::memcpy(datap + datalen_seqlen, & length, sizeof(unsigned int));

      if (length == 0)
        {
          fatal(error_prefix, "Empty sequence found on line ", lineno - 1, ".");
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

          std::memcpy(datap + datalen, & nt_buffer, sizeof(nt_buffer));
          datalen += sizeof(nt_buffer);

          nt_buffer = 0;
          nt_bufferlen = 0;  // that value is never read again, all tests pass without it
        }

      ++sequences;

      if (is_regular) {
        progress_update(filepos);
      }
    }
  progress_done();

  std::fclose(input_fp);

  /* init zobrist hashing */

  // add 2 for two insertions (refactoring: insertions in headers?)
  const unsigned int zobrist_len = std::max(4 * longestheader, longest + 2);
  zobrist_init(zobrist_len);

  /* set up hash to check for unique headers */

  const uint64_t hdrhashsize {2 * sequences};

  auto * * hdrhashtable = new seqinfo_t*[hdrhashsize] { };

  /* set up hash to check for unique sequences */

  const uint64_t seqhashsize {2 * sequences};

  seqinfo_t * * seqhashtable {nullptr};

  if (parameters.opt_differences > 1)
    {
      seqhashtable = new seqinfo_t*[seqhashsize] { };
    }

  /* create indices */

  seqindex = new seqinfo_t[sequences];
  seqinfo_t * seqindex_p {seqindex};

  seqinfo_t * lastseq {nullptr};

  bool presorted {true};

  char * pl {datap};  // refactoring: purpose and how to rename?
  progress_init("Indexing database:", sequences);
  for(auto i = 0ULL; i < sequences; i++)
    {
      /* get line number */
      const unsigned int line_number = *(reinterpret_cast<unsigned int*>(pl));
      pl += sizeof(unsigned int);

      /* get header */
      seqindex_p->header = pl;
      seqindex_p->headerlen = static_cast<int>(std::strlen(seqindex_p->header));
      pl += seqindex_p->headerlen + 1;

      /* and sequence */
      const unsigned int seqlen = *(reinterpret_cast<unsigned int*>(pl));
      seqindex_p->seqlen = seqlen;
      pl += sizeof(unsigned int);
      seqindex_p->seq = pl;
      pl += nt_bytelength(seqlen);

      /* get amplicon abundance */
      find_abundance(seqindex_p, line_number, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      if ((seqindex_p->abundance_start == 0) &&
          (seqindex_p->abundance_end == seqindex_p->headerlen)) {
        fatal(error_prefix, "Empty sequence identifier.");
      }

      /* check if the sequences are presorted by abundance and header */

      if (presorted && (lastseq != nullptr))
        {
          if (lastseq->abundance < seqindex_p->abundance) {
            presorted = false;
          }
          else if (lastseq->abundance == seqindex_p->abundance)
            {
              if (std::strcmp(lastseq->header, seqindex_p->header) > 0) {
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

      const uint64_t hdrhash = zobrist_hash(reinterpret_cast<unsigned char*>
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
                  (std::strncmp(seqindex_p->header + id_start,
                                hdrfound->header + hit_id_start,
                                static_cast<uint64_t>(id_len)) == 0)) {
                break;
              }
            }

          hdrhashindex = (hdrhashindex + 1) % hdrhashsize;
        }

      if (hdrfound != nullptr)
        {
          const std::string full_header {seqindex_p->header + id_start};
          fatal(error_prefix, "Duplicated sequence identifier: ", full_header.substr(0, id_len));
        }

      hdrhashtable[hdrhashindex] = seqindex_p;

      /* hash sequence */
      seqindex_p->seqhash = zobrist_hash(reinterpret_cast<unsigned char*>
                                         (seqindex_p->seq),
                                         seqindex_p->seqlen);

      if (parameters.opt_differences > 1)
        {
          /* Check for duplicated sequences using hash table, */
          /* but only for d>1. Handled internally for d=1.    */

          uint64_t seqhashindex = seqindex_p->seqhash % seqhashsize;
          seqinfo_t * seqfound {nullptr};

          while ((seqfound = seqhashtable[seqhashindex]) != nullptr)
            {
              if ((seqfound->seqhash == seqindex_p->seqhash) &&
                  (seqfound->seqlen == seqindex_p->seqlen) &&
                  (std::memcmp(seqfound->seq,
                          seqindex_p->seq,
                          nt_bytelength(seqindex_p->seqlen)) == 0)) {
                break;
              }
              seqhashindex = (seqhashindex + 1) % seqhashsize;
            }

          if (seqfound != nullptr)
            {
              ++duplicates_found;
              break;
            }
          seqhashtable[seqhashindex] = seqindex_p;
        }

      ++seqindex_p;
      progress_update(i);
    }

  if (duplicates_found != 0U)
    {
      fatal(error_prefix,
            "some fasta entries have identical sequences.\n"
            "Swarm expects dereplicated fasta files.\n"
            "Such files can be produced with swarm or vsearch:\n"
            " swarm -d 0 -w derep.fasta -o /dev/null input.fasta\n"
            "or\n"
            " vsearch --derep_fulllength input.fasta --sizein --sizeout --output derep.fasta");
    }

  progress_done();

  delete [] line;
  line = nullptr;
  linecap = 0;

  if (missingabundance != 0)
    {
      fatal(error_prefix, "Abundance annotations not found for ",
            missingabundance, " sequences, starting on line ",
            missingabundance_lineno, ".\n>",
            missingabundance_header, "\n",
            "Fasta headers must end with abundance annotations (_INT or ;size=INT).\n"
            "The -z option must be used if the abundance annotation is in the latter format.\n"
            "Abundance annotations can be produced by dereplicating the sequences.\n"
            "The header is defined as the string comprised between the \">\" symbol\n"
            "and the first space or the end of the line, whichever comes first.");
    }

  if (! presorted)
    {
      progress_init("Abundance sorting:", 1);
      std::qsort(seqindex, sequences, sizeof(seqinfo_t), db_compare_abundance);
      progress_done();
    }

  delete [] hdrhashtable;
  hdrhashtable = nullptr;

  delete [] seqhashtable;
  seqhashtable = nullptr;

  // user report
  std::fprintf(logfile, "Database info:     %" PRIu64 " nt", db_getnucleotidecount());
  std::fprintf(logfile, " in %u sequences,", db_getsequencecount());
  std::fprintf(logfile, " longest %u nt\n", db_getlongestsequence());
}


void db_qgrams_init()
{
  qgrams = new qgramvector_t[sequences];

  seqinfo_t * seqindex_p {seqindex};
  progress_init("Find qgram vects: ", sequences);
  for(auto i = 0U; i < sequences; i++)
    {
      /* find qgrams */
      findqgrams(reinterpret_cast<unsigned char*>(seqindex_p->seq),
                 seqindex_p->seqlen,
                 qgrams[i]);
      ++seqindex_p;
      progress_update(i);
    }
  progress_done();
}


void db_qgrams_done()
{
  delete [] qgrams;
  qgrams = nullptr;
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
  datap = nullptr;
  delete [] seqindex;
  seqindex = nullptr;
}


auto db_fprintseq(std::FILE * fastaout_fp, const unsigned int seqno) -> void
{
  const unsigned int len {db_getsequencelen(seqno)};
  char * const seqptr {db_getsequence(seqno)};
  static std::vector<char> buffer(db_getlongestsequence() + 1, '\0');

  // decode to nucleotides (A, C, G and T)
  for(auto i = 0U; i < len; i++) {
    buffer[i] = sym_nt[1 + nt_extract(seqptr, i)];
  }
  buffer[len] = '\0';

  std::fprintf(fastaout_fp, "%.*s\n", len, buffer.data());
}
