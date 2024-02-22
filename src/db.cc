/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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
#include "utils/input_output.h"
#include "utils/nt_codec.h"
#include "utils/progress.h"
#include "utils/qgram_array.h"
#include "utils/seqinfo.h"
#include "zobrist.h"
#include <algorithm>  // std::max() std::min() std::sort()
#include <array>
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fileno, fclose(), size_t // stdio.h: fdopen, ssize_t, getline
#include <cstdlib>  // qsort()
#include <cstring>  // memcpy, memcmp
#include <iterator>  // std::next()
#include <limits>
#include <string>
#include <sys/stat.h>  // fstat, S_ISREG, stat
#include <vector>

#ifndef PRIu64
#ifdef _WIN32
#define PRIu64 "I64u"
#else
constexpr char PRIu64[] = "lu";
#endif
#endif


constexpr unsigned int memchunk {1U << 20U};  // 1 megabyte
constexpr unsigned int linealloc {2048};
constexpr auto int8_max = std::numeric_limits<int8_t>::max();
constexpr long unsigned int n_chars {int8_max + 1};  // 128 ascii chars
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


struct Seq_stats {
  uint64_t nucleotides {0};
  unsigned int longestheader {0};
  int missingabundance {0};
  uint64_t missingabundance_lineno {0};
  char * missingabundance_header {nullptr};
};

static unsigned int sequences {0};
static unsigned int longest {0};

struct seqinfo_s * seqindex {nullptr};


auto db_getsequencecount() -> unsigned int
{
  return sequences;
}


auto db_getlongestsequence() -> unsigned int
{
  return longest;
}


auto fprint_id(std::FILE * stream, const uint64_t seqno, const bool opt_usearch_abundance,
               const int64_t opt_append_abundance) -> void
{
  auto const & seqinfo = *std::next(seqindex, static_cast<std::ptrdiff_t>(seqno));
  auto const * hdrstr = seqinfo.header;
  auto const hdrlen = seqinfo.headerlen;
  auto const abundance = seqinfo.abundance;

  // if abundance is missing and if user says that a missing abundance is ok, then...
  if ((opt_append_abundance != 0) and (seqinfo.abundance_start == seqinfo.abundance_end)) {
    if (opt_usearch_abundance) {
      std::fprintf(stream, "%.*s;size=%" PRIu64 ";", hdrlen, hdrstr, abundance);
    }
    else {
      std::fprintf(stream, "%.*s_%" PRIu64, hdrlen, hdrstr, abundance);
    }
  }
  else {
    std::fprintf(stream, "%.*s", hdrlen, hdrstr);
  }
}


auto fprint_id_noabundance(std::FILE * stream, const uint64_t seqno, const bool opt_usearch_abundance) -> void
{
  auto const & seqinfo = *std::next(seqindex, static_cast<std::ptrdiff_t>(seqno));
  auto const * hdrstr = seqinfo.header;
  auto const hdrlen = seqinfo.headerlen;
  auto const abundance_start = seqinfo.abundance_start;
  auto const abundance_end = seqinfo.abundance_end;

  if (abundance_start < abundance_end)
    {
      /* print start of header */
      std::fprintf(stream, "%.*s", abundance_start, hdrstr);

      if (opt_usearch_abundance)
        {
          /* print semicolon if the abundance is not at either end */
          if ((abundance_start > 0) and (abundance_end < hdrlen)) {
            std::fprintf(stream, ";");
          }

          /* print remaining part */
          std::fprintf(stream, "%.*s", hdrlen - abundance_end, std::next(hdrstr, abundance_end));
        }
    }
  else {
    std::fprintf(stream, "%.*s", hdrlen, hdrstr);
  }
}


auto fprint_id_with_new_abundance(std::FILE * stream,
                                  const uint64_t seqno,
                                  const uint64_t abundance,
                                  const bool opt_usearch_abundance) -> void
{
  auto const & seqinfo = *std::next(seqindex, static_cast<std::ptrdiff_t>(seqno));

  if (opt_usearch_abundance) {
    std::fprintf(stream,
                 "%.*s%ssize=%" PRIu64 ";%.*s",
                 seqinfo.abundance_start,
                 seqinfo.header,
                 seqinfo.abundance_start > 0 ? ";" : "",
                 abundance,
                 seqinfo.headerlen - seqinfo.abundance_end,
                 std::next(seqinfo.header, seqinfo.abundance_end));
  }
  else {
    std::fprintf(stream,
                 "%.*s_%" PRIu64,
                 seqinfo.abundance_start,
                 seqinfo.header,
                 abundance);
  }
}


auto find_swarm_abundance(const char * header,
                          int & start,
                          int & end,
                          int64_t & number) -> bool
{
  /*
    Identify the first occurence of the pattern (_)([0-9]+)$
    in the header string.
  */

  start = 0;
  end = 0;
  number = 0;

  static constexpr unsigned int max_digits {20};  // 20 digits at most (abundance > 10^20)
  static const std::string digit_chars = "0123456789";

  assert(header != nullptr); // assert to prove impossible
  if (header == nullptr) {
    return false;  // refactoring: if header cannot be a nullptr, replace with assert
  }

  const char * abundance_string = std::strrchr(header, '_');

  if (abundance_string == nullptr) {
    return false;
  }

  const std::size_t n_digits = std::strspn(std::next(abundance_string), digit_chars.c_str());

  if (n_digits > max_digits) {
    return false;
  }

  if (*std::next(abundance_string, static_cast<std::ptrdiff_t>(n_digits + 1)) != 0) {
    return false;
  }

  const int64_t abundance_start = abundance_string - header;
  const int64_t abundance_end = abundance_start + 1 + static_cast<int64_t>(n_digits);

  start = static_cast<int>(abundance_start);
  end = static_cast<int>(abundance_end);
  number = std::atol(std::next(abundance_string)); // refactoring: std::strtol(start, end, base)

  return true;
}


auto find_usearch_abundance(const char * header,
                            int & start,
                            int & end,
                            int64_t & number) -> bool
{
  /*
    Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
    in the header string.
  */

  assert(header != nullptr); // header cannot be a nullptr at this stage

  static const std::string attribute {"size="};
  static const std::string digit_chars {"0123456789"};
  auto const hlen = static_cast<int64_t>(std::strlen(header));
  auto const alen = static_cast<int64_t>(attribute.length());
  int64_t position = 0;

  while (position + alen < hlen)
    {
      auto const * result = std::strstr(std::next(header, position), attribute.c_str());

      /* no match */
      assert(result != nullptr); // assert to prove impossible
      if (result == nullptr) {
        break;
      }

      position = result - header;

      /* check for ';' in front */
      if ((position > 0) and (*std::next(header, position - 1) != ';'))
        {
          position += alen + 1;
          continue;
        }

      auto const n_digits = static_cast<int64_t>(std::strspn(std::next(header, position + alen), digit_chars.c_str()));

      /* check for at least one digit */
      if (n_digits == 0)
        {
          position += alen + 1;
          continue;
        }

      /* check for ';' after */
      if ((position + alen + n_digits < hlen) and (*std::next(header, position + alen + n_digits) != ';'))
        {
          position += alen + n_digits + 2;
          continue;
        }

      /* ok */
      if (position > 0) {
        start = static_cast<int>(position - 1);
      }
      else {
        start = 0;
      }
      end   = static_cast<int>(std::min(position + alen + n_digits + 1, hlen));
      number = std::atol(std::next(header, position + alen));

      return true;
    }

  return false;
}


auto find_abundance(struct seqinfo_s & seqinfo, struct Seq_stats & seq_stats, uint64_t lineno,
                    bool opt_usearch_abundance, int64_t opt_append_abundance) -> void
{
  char * header = seqinfo.header;

  /* read size/abundance annotation */
  int64_t abundance = 0;
  int start = 0;
  int end = 0;
  int64_t number = 0;

  if (opt_usearch_abundance)
    {
      /* (^|;)size=([0-9]+)(;|$) */

      if (find_usearch_abundance(header, start, end, number))
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

      if (find_swarm_abundance(header, start, end, number))
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
      start = seqinfo.headerlen;
      end = start;

      if (opt_append_abundance != 0) {
        abundance = opt_append_abundance;
      }
      else
        {
          ++seq_stats.missingabundance;
          // record the position of the first missing abundance entry
          if (seq_stats.missingabundance == 1)
            {
              seq_stats.missingabundance_lineno = lineno;
              seq_stats.missingabundance_header = header;
            }
        }
    }

  seqinfo.abundance = static_cast<uint64_t>(abundance);
  seqinfo.abundance_start = start;
  seqinfo.abundance_end = end;
}


auto sort_index_if_need_be(struct Parameters const & parameters,
                           std::vector<struct seqinfo_s> & seqindex_v) -> void {
      progress_init("Abundance sorting:", 1);

      auto compare_entries = [](struct seqinfo_s const& lhs,
                                struct seqinfo_s const& rhs) -> bool
      {
        // sort by decreasing abundance
        if (lhs.abundance > rhs.abundance) {
          return true;
        }

        if (lhs.abundance < rhs.abundance) {
          return false;
        }

        // ...then ties are sorted by header (lexicographical order)
        return std::strcmp(lhs.header, rhs.header) < 0;
      };

      if (not std::is_sorted(seqindex_v.begin(), seqindex_v.end(),
                             compare_entries)) {
        std::sort(seqindex_v.begin(), seqindex_v.end(), compare_entries);
      }
      progress_done(parameters);
}


auto db_read(struct Parameters const & parameters,
             std::vector<char> & data_v,
             std::vector<struct seqinfo_s> & seqindex_v,
             std::vector<uint64_t> & zobrist_tab_base_v,
             std::vector<uint64_t> & zobrist_tab_byte_base_v) -> void
{
  struct Seq_stats seq_stats;
  uint64_t datalen {0};
  uint64_t duplicates_found {0};

  longest = 0;
  sequences = 0;

  /* open input file or stream */

  assert(parameters.input_filename.c_str() != nullptr);  // filename is set to '-' (stdin) by default

  std::FILE * input_fp { fopen_input(parameters.input_filename.c_str()) };
  if (input_fp == nullptr)
    {
      fatal(error_prefix, "Unable to open input data file (", parameters.input_filename.c_str(), ").\n");
    }

  /* get file size */
  // refactoring: C++17 std::filesystem::file_size
  struct stat fstat_buffer;  // refactoring: add initializer '{}' (warning with GCC < 5)

  if (fstat(fileno(input_fp), &fstat_buffer) != 0)  // refactor: fstat and fileno linuxisms
    {
      fatal(error_prefix, "Unable to fstat on input file (", parameters.input_filename.c_str(), ").\n");
    }
  const bool is_regular = S_ISREG(fstat_buffer.st_mode);  // refactoring: S_ISREG linuxisms
  const uint64_t filesize = is_regular ? static_cast<uint64_t>(fstat_buffer.st_size) : 0;
  uint64_t filepos = 0;

  if (not is_regular) {
    std::fprintf(parameters.logfile, "Waiting for data... (hit Ctrl-C and run 'swarm -h' if you meant to read data from a file)\n");
  }

  std::size_t linecap = linealloc;
  auto * line = static_cast<char *>(xmalloc(linecap)); // char * line {new char[linecap]};  // refactoring: replacing with a std::vector fails, as getline might need to reallocate and will free() 'line', creating a double-free attempt at the end of the scope
  ssize_t linelen = xgetline(& line, & linecap, input_fp);
  if (linelen < 0)
    {
      line[0] = 0;
      linelen = 0;
    }
  filepos += static_cast<unsigned long int>(linelen);

  auto lineno = 1U;

  progress_init("Reading sequences:", filesize);

  /* allocate space */
  if (filesize > memchunk) {
    // in-RAM data cannot be smaller than 1/4 of the on-disk data
    data_v.reserve(filesize / 4);
  }
  data_v.resize(memchunk);

  while(line[0] != 0)
    {
      /* read header */
      /* the header ends at a space, cr, lf or null character */

      if (line[0] != '>') {
        fatal(error_prefix, "Illegal header line in fasta file.");
      }

      auto headerlen = static_cast<unsigned int>
        (std::strcspn(std::next(line), " \r\n"));

      seq_stats.longestheader = std::max(headerlen, seq_stats.longestheader);

      if (seq_stats.longestheader > max_header_length) {
        fatal(error_prefix, "Headers longer than 16,777,215 symbols are not supported.");
      }

      /* store the line number */

      while (datalen + sizeof(unsigned int) > data_v.size())
        {
          data_v.resize(data_v.size() + memchunk);
        }
      std::memcpy(&data_v[datalen], & lineno, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* store the header */

      while (datalen + headerlen + 1 > data_v.size())
        {
          data_v.resize(data_v.size() + memchunk);
        }
      std::memcpy(&data_v[datalen], std::next(line), headerlen);
      data_v[datalen + headerlen] = 0;
      datalen += headerlen + 1;


      /* get next line */

      linelen = xgetline(& line, & linecap, input_fp);
      if (linelen < 0)
        {
          line[0] = 0;
          linelen = 0;
        }
      filepos += static_cast<unsigned long int>(linelen);

      ++lineno;


      /* store a dummy sequence length */

      auto length = 0U;

      while (datalen + sizeof(unsigned int) > data_v.size())
        {
          data_v.resize(data_v.size() + memchunk);
        }
      const uint64_t datalen_seqlen = datalen;
      std::memcpy(&data_v[datalen], & length, sizeof(unsigned int));
      datalen += sizeof(unsigned int);


      /* read and store sequence */

      uint64_t nt_buffer {0};
      auto nt_bufferlen = 0U;
      static constexpr unsigned int nt_buffersize {4 * sizeof(nt_buffer)};
      static constexpr unsigned char null_char = '\0';
      static constexpr int new_line {10};
      static constexpr int carriage_return {13};
      static constexpr int start_chars_range {32};  // visible ascii chars: 32-126
      static constexpr int end_chars_range {126};

      while ((line[0] != 0) and (line[0] != '>'))
        {
          auto character = null_char;
          auto * line_ptr = line;
          while((character = static_cast<unsigned char>(*line_ptr)) != null_char)
            {
              line_ptr = std::next(line_ptr);
              const auto mapped_char = map_nt[character];
              if (mapped_char != 0)
                {
                  nt_buffer |= (mapped_char - 1) << (2 * nt_bufferlen);
                  ++length;
                  ++nt_bufferlen;

                  if (nt_bufferlen == nt_buffersize)
                    {
                      while (datalen + sizeof(nt_buffer) > data_v.size())
                        {
                          data_v.resize(data_v.size() + memchunk);
                        }

                      std::memcpy(&data_v[datalen], & nt_buffer, sizeof(nt_buffer));
                      datalen += sizeof(nt_buffer);

                      nt_bufferlen = 0;
                      nt_buffer = 0;
                    }
                }
              else if ((character != new_line) and (character != carriage_return))
                {
                  if ((character >= start_chars_range) and (character <= end_chars_range)) {
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
          filepos += static_cast<unsigned long int>(linelen);

          ++lineno;
        }

      /* fill in real length */

      std::memcpy(&data_v[datalen_seqlen], & length, sizeof(unsigned int));

      if (length == 0)
        {
          fatal(error_prefix, "Empty sequence found on line ", lineno - 1, ".");
        }

      seq_stats.nucleotides += length;
      longest = std::max(length, longest);


      /* save remaining padded 64-bit value with nt's, if any */

      if (nt_bufferlen > 0)
        {
          while (datalen + sizeof(nt_buffer) > data_v.size())
            {
              data_v.resize(data_v.size() + memchunk);
            }

          std::memcpy(&data_v[datalen], & nt_buffer, sizeof(nt_buffer));
          datalen += sizeof(nt_buffer);

          nt_buffer = 0;
          nt_bufferlen = 0;  // that value is never read again, all tests pass without it
        }

      ++sequences;

      if (is_regular) {
        progress_update(filepos);
      }
    }
  progress_done(parameters);

  std::fclose(input_fp);

  /* init zobrist hashing */

  // add 2 for two insertions (refactoring: insertions in headers?)
  const auto zobrist_len = std::max(4 * seq_stats.longestheader, longest + 2);
  zobrist_init(zobrist_len, zobrist_tab_base_v, zobrist_tab_byte_base_v);

  /* set up hash to check for unique headers */

  const uint64_t hdrhashsize {2ULL * sequences};
  std::vector<struct seqinfo_s *> hdrhashtable(hdrhashsize);

  /* set up hash to check for unique sequences */

  const uint64_t seqhashsize {2ULL * sequences};

  std::vector<struct seqinfo_s *> seqhashtable;

  if (parameters.opt_differences > 1) {
    seqhashtable.resize(seqhashsize);
  }

  /* create indices */

  seqindex_v.resize(sequences);
  seqindex = seqindex_v.data();

  auto * cursor = data_v.data();
  progress_init("Indexing database:", sequences);
  auto counter = 0ULL;
  for(auto& a_sequence: seqindex_v) {

      /* get line number */
      const auto line_number = *(reinterpret_cast<unsigned int*>(cursor));  // UBSAN: misaligned address for type 'unsigned int', which requires 4 byte alignment
      cursor = std::next(cursor, sizeof(unsigned int));

      /* get header */
      a_sequence.header = cursor;
      a_sequence.headerlen = static_cast<int>(std::strlen(a_sequence.header));
      cursor = std::next(cursor, a_sequence.headerlen + 1);

      /* and sequence */
      const auto seqlen = *(reinterpret_cast<unsigned int*>(cursor));  // UBSAN: misaligned address for type 'unsigned int', which requires 4 byte alignment
      a_sequence.seqlen = seqlen;
      cursor = std::next(cursor, sizeof(unsigned int));
      a_sequence.seq = cursor;
      cursor = std::next(cursor, nt_bytelength(seqlen));

      /* get amplicon abundance */
      find_abundance(a_sequence, seq_stats, line_number, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

      if ((a_sequence.abundance_start == 0) and
          (a_sequence.abundance_end == a_sequence.headerlen)) {
        fatal(error_prefix, "Empty sequence identifier.");
      }

      /* check for duplicated identifiers using hash table */
      // refactoring: extract to a free function, perform for each new header
      // C++14 refactoring: std::set::find() heterogeneous lookup (see overloads 3 and 4,
      // https://en.cppreference.com/w/cpp/container/set/find)

      /* find position and length of identifier in header */

      int id_start {0};
      int id_len {0};

      if (a_sequence.abundance_start > 0)
        {
          /* id first, then abundance (e.g. >name;size=1 or >name_1) */
          id_start = 0;
          id_len = a_sequence.abundance_start;
        }
      else
        {
          /* abundance first then id (e.g. >size=1;name) */
          id_start = a_sequence.abundance_end;
          id_len = a_sequence.headerlen - a_sequence.abundance_end;
        }

      const auto hdrhash = zobrist_hash(reinterpret_cast<unsigned char*>
                                        (std::next(a_sequence.header, id_start)),
                                        4 * static_cast<unsigned int>(id_len));

      a_sequence.hdrhash = hdrhash;
      uint64_t hdrhashindex = hdrhash % hdrhashsize;

      struct seqinfo_s * hdrfound {nullptr};

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

              if ((id_len == hit_id_len) and
                  (std::strncmp(std::next(a_sequence.header, id_start),
                                std::next(hdrfound->header, hit_id_start),
                                static_cast<uint64_t>(id_len)) == 0)) {
                break;
              }
            }

          hdrhashindex = (hdrhashindex + 1) % hdrhashsize;
        }

      if (hdrfound != nullptr)
        {
          const std::string full_header {std::next(a_sequence.header, id_start)};
          fatal(error_prefix, "Duplicated sequence identifier: ", full_header.substr(0, static_cast<unsigned long int>(id_len)));
        }

      hdrhashtable[hdrhashindex] = &a_sequence;

      /* hash sequence */
      a_sequence.seqhash = zobrist_hash(reinterpret_cast<unsigned char*>
                                         (a_sequence.seq),
                                         a_sequence.seqlen);

      if (parameters.opt_differences > 1)
        {
          // refactoring: extract to a free function (not trivial)
          /* Check for duplicated sequences using hash table,  */
          /* but only for d > 1. Handled internally for d = 1. */

          uint64_t seqhashindex = a_sequence.seqhash % seqhashsize;
          struct seqinfo_s * seqfound {nullptr};

          while ((seqfound = seqhashtable[seqhashindex]) != nullptr)
            {
              if ((seqfound->seqhash == a_sequence.seqhash) and
                  (seqfound->seqlen == a_sequence.seqlen) and
                  std::equal(seqfound->seq,
                             std::next(seqfound->seq, nt_bytelength(a_sequence.seqlen)),
                             a_sequence.seq)) {
                break;
              }
              seqhashindex = (seqhashindex + 1) % seqhashsize;
            }

          if (seqfound != nullptr)
            {
              ++duplicates_found;
              break;
            }
          seqhashtable[seqhashindex] = &a_sequence;
        }

      progress_update(counter);
      ++counter;
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

  progress_done(parameters);

  if (line != nullptr)
    {
      xfree(line);
      line = nullptr;
      linecap = 0;
    }

  if (seq_stats.missingabundance != 0)
    {
      fatal(error_prefix, "Abundance annotations not found for ",
            seq_stats.missingabundance, " sequences, starting on line ",
            seq_stats.missingabundance_lineno, ".\n>",
            seq_stats.missingabundance_header, "\n",
            "Fasta headers must end with abundance annotations (_INT or ;size=INT).\n"
            "The -z option must be used if the abundance annotation is in the latter format.\n"
            "Abundance annotations can be produced by dereplicating the sequences.\n"
            "The header is defined as the string comprised between the \">\" symbol\n"
            "and the first space or the end of the line, whichever comes first.");
    }

  sort_index_if_need_be(parameters, seqindex_v);

  // user report
  std::fprintf(parameters.logfile, "Database info:     %" PRIu64 " nt", seq_stats.nucleotides);
  std::fprintf(parameters.logfile, " in %u sequences,", db_getsequencecount());
  std::fprintf(parameters.logfile, " longest %u nt\n", db_getlongestsequence());
}


// refactoring: only used in algo.cc, extract to its own header file?
auto db_qgrams_init(struct Parameters const & parameters,
                    std::vector<struct seqinfo_s> & seqindex_v) -> void
{
  // refactoring: qgrams is an array of char arrays!
  // - vector of std::array is not allowed,
  // - vector of vector is not contiguous!
  // - is contiguity really a requirement??
  // in the meantime:
  // - std::vector<char> qgrams_v(unitSize * sequences, '\0');  // unitSize = qgramvectorbytes = 128
  // - or std::vector<std::vector<char>> qgrams_v(sequences, std::vector<char>(unitSize, '\0'));
  qgrams = new qgramvector_t[sequences];

  progress_init("Find qgram vects: ", seqindex_v.size());
  auto counter = 0U;
  for(auto const & seqindex_p : seqindex_v) {
    /* find qgrams */
    findqgrams(reinterpret_cast<unsigned char*>(seqindex_p.seq),
               seqindex_p.seqlen,
               qgrams[counter]);
    progress_update(counter);
    ++counter;
  }
  progress_done(parameters);
}


auto db_qgrams_done() -> void
{
  delete [] qgrams;
  qgrams = nullptr;
}


auto db_gethash(const uint64_t seqno) -> uint64_t
{
  auto const record_number = static_cast<std::ptrdiff_t>(seqno);
  auto const & fasta_record = *std::next(seqindex, record_number);
  return fasta_record.seqhash;
}


auto db_getsequence(const uint64_t seqno) -> char *
{
  auto const record_number = static_cast<std::ptrdiff_t>(seqno);
  auto const & fasta_record = *std::next(seqindex, record_number);
  return fasta_record.seq;
}


auto db_getsequenceandlength(uint64_t seqno,
                             char * & address,
                             unsigned int & length) -> void
{
  auto const record_number = static_cast<std::ptrdiff_t>(seqno);
  auto const & fasta_record = *std::next(seqindex, record_number);
  address = fasta_record.seq;
  length = fasta_record.seqlen;
}


auto db_getsequencelen(const uint64_t seqno) -> unsigned int
{
  auto const record_number = static_cast<std::ptrdiff_t>(seqno);
  auto const & fasta_record = *std::next(seqindex, record_number);
  return fasta_record.seqlen;
}


auto db_getheader(const uint64_t seqno) -> char *
{
  auto const record_number = static_cast<std::ptrdiff_t>(seqno);
  auto const & fasta_record = *std::next(seqindex, record_number);
  return fasta_record.header;
}


auto db_getabundance(const uint64_t seqno) -> uint64_t
{
  auto const record_number = static_cast<std::ptrdiff_t>(seqno);
  auto const & fasta_record = *std::next(seqindex, record_number);
  return fasta_record.abundance;
}


auto db_free() -> void
{
  seqindex = nullptr;
}


auto db_fprintseq(std::FILE * fastaout_fp, const unsigned int seqno) -> void
{
  auto const len = db_getsequencelen(seqno);
  auto * const seqptr = db_getsequence(seqno);
  static std::vector<char> buffer(db_getlongestsequence() + 1, '\0');

  // decode to nucleotides (A, C, G and T)
  for(auto i = 0U; i < len; i++) {
    buffer[i] = sym_nt[1 + nt_extract(seqptr, i)];
  }
  buffer[len] = '\0';

  std::fprintf(fastaout_fp, "%.*s\n", len, buffer.data());
}
