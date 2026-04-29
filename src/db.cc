/*
    SWARM

    Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

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
#include "db.h"
#include "utils/fatal.h"
#include "utils/input_output.h"
#include "utils/nt_codec.h"
#include "utils/progress.h"
#include "utils/seq_index.h"
#include "utils/view.h"
#include "utils/xgetline.h"
#include <algorithm>  // std::max() std::min() std::sort()
#include <array>
#include <cassert>  // assert()
#include <cinttypes>  // macros PRIu64 and PRId64
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // fileno, size_t // stdio.h: fdopen, ssize_t, getline
#include <cstdlib>  // qsort()
#include <cstring>  // memcpy
#include <iterator>  // std::next()
#include <limits>
#include <memory>  // std::unique_ptr
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


// anonymous namespace: limit visibility and usage to this translation unit
namespace {

  constexpr unsigned int memchunk {1U << 20U};  // 1 megabyte
  constexpr auto int8_max = std::numeric_limits<int8_t>::max();
  constexpr long unsigned int n_chars {int8_max + 1};  // 128 ascii chars

  struct File_info {
    uint64_t filesize {0};
    bool is_regular {false};
  };


  struct Seq_stats {
    uint64_t nucleotides {0};
    unsigned int longestheader {0};
    int missingabundance {0};
    uint64_t missingabundance_lineno {0};
    char const * missingabundance_header {nullptr};
    unsigned int n_sequences {0};
    unsigned int longest_sequence {0};
    bool has_duplicates {false};
  };


  struct Parse_result {
    std::vector<struct Entry> entries;
    struct Seq_stats stats;
  };


  // RAII wrapper for the line buffer passed to xgetline(). POSIX
  // getline() owns the buffer's lifetime: it may std::realloc() it on
  // long lines, so the storage must come from std::malloc and the
  // destructor must call std::free. std::vector<char> or new[]/delete[]
  // would create an allocator mismatch and undefined behavior.
  struct Line_buffer {
    char *      data {nullptr};
    std::size_t capacity {0};

    explicit Line_buffer(std::size_t const initial)
      : data{static_cast<char *>(std::malloc(initial))}, capacity{initial}
    {
      if (data == nullptr) {
        fatal(error_prefix, "Unable to allocate enough memory.");
      }
    }

    // noexcept: std::free is noexcept and the nullptr guard performs
    // only an integer comparison.
    ~Line_buffer() noexcept { release(); }

    // noexcept: see destructor.
    auto release() noexcept -> void {
      if (data != nullptr) {
        std::free(data);
        data = nullptr;
        capacity = 0;
      }
    }

    Line_buffer(Line_buffer const &)                     = delete;
    auto operator=(Line_buffer const &) -> Line_buffer & = delete;
    Line_buffer(Line_buffer &&)                          = delete;
    auto operator=(Line_buffer &&)      -> Line_buffer & = delete;
  };


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


  auto warn_if_file_is_not_regular(struct Parameters const & parameters, bool const is_regular) -> void {
    if (not is_regular) {
      std::fprintf(parameters.logfile, "Waiting for data... (hit Ctrl-C and run 'swarm -h' if you meant to read data from a file)\n");
    }
  }


  auto get_file_info(std::FILE * input_handle, struct Parameters const & parameters) -> struct File_info {
    // get file size and file type (regular or pipe)
    // refactoring: C++17 std::filesystem::file_size
    struct File_info file_info;;
    struct stat fstat_buffer;  // refactoring: add initializer '{}' (warning with GCC < 5)

    if (fstat(fileno(input_handle), &fstat_buffer) != 0) { // refactor: fstat and fileno are linuxisms
      fatal(error_prefix, "Unable to fstat on input file (", parameters.input_filename.c_str(), ").\n");
    }
    file_info.is_regular = S_ISREG(fstat_buffer.st_mode);  // refactoring: S_ISREG is a linuxism
    file_info.filesize = file_info.is_regular ? static_cast<uint64_t>(fstat_buffer.st_size) : 0U;
    return file_info;
  }


  auto initial_allocation(std::vector<char> & data_v, uint64_t const filesize) -> void {
    auto const minimal_reserve = filesize >> 2U;  // 1/4 of filesize
    if (minimal_reserve > memchunk) {
      // in-RAM data cannot be smaller than 1/4 of the on-disk data
      data_v.reserve(minimal_reserve);
    }
    data_v.resize(memchunk);
  }


  auto linear_resize_if_need_be(std::vector<char> & data_v, uint64_t const minimal_size) -> void {
    auto const current_size = data_v.size();
    if (current_size > minimal_size) { return; }
    auto new_size = current_size;
    while (minimal_size > new_size) {
      assert(new_size <= std::numeric_limits<uint64_t>::max() - memchunk);
      new_size += memchunk;
    }
    data_v.resize(new_size);
  }


  auto find_swarm_abundance(View<char const> const header_view,
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

    // strrchr / strspn / strtoll require a null-terminated string;
    // header_view always points into Data::data_, where each header
    // is followed by a '\0' byte written at parse time.
    auto const * const header = header_view.data();

    assert(header != nullptr); // assert to prove impossible
    if (header == nullptr) {
      return false;  // refactoring: if header cannot be a nullptr, replace with assert
    }

    auto const * const abundance_string = std::strrchr(header, '_');

    if (abundance_string == nullptr) {
      return false;
    }

    std::size_t const n_digits = std::strspn(std::next(abundance_string), digit_chars.c_str());

    if (n_digits > max_digits) {
      return false;
    }

    assert((n_digits + 1) <= std::numeric_limits<std::ptrdiff_t>::max());
    if (*std::next(abundance_string, static_cast<std::ptrdiff_t>(n_digits + 1)) != 0) {
      return false;
    }

    int64_t const abundance_start = std::distance(header_view.cbegin(), abundance_string);
    assert(n_digits <= std::numeric_limits<int64_t>::max());
    int64_t const abundance_end = abundance_start + 1 + static_cast<int64_t>(n_digits);

    assert(abundance_start <= std::numeric_limits<int>::max());
    assert(abundance_end <= std::numeric_limits<int>::max());
    start = static_cast<int>(abundance_start);
    end = static_cast<int>(abundance_end);
    // refactoring: capture strtoll's end pointer and check errno == ERANGE
    // to detect overflow (n_digits is bounded above by max_digits = 20,
    // which can exceed int64_t's 19-digit range).
    static constexpr int base_value {10};
    number = std::strtoll(std::next(abundance_string), nullptr, base_value);

    return true;
  }


  auto find_usearch_abundance(View<char const> const header_view,
                              int & start,
                              int & end,
                              int64_t & number) -> bool
  {
    /*
      Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
      in the header string.
    */

    // strstr / strspn / strtoll require a null-terminated string;
    // header_view always points into Data::data_, where each header
    // is followed by a '\0' byte written at parse time.
    auto const * const header = header_view.data();
    assert(header != nullptr); // header cannot be a nullptr at this stage

    static const std::string attribute {"size="};
    static const std::string digit_chars {"0123456789"};
    auto const hlen = static_cast<int64_t>(header_view.size());
    assert(attribute.length() <= std::numeric_limits<int64_t>::max());
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
          assert((position - 1) <= std::numeric_limits<int>::max());
          start = static_cast<int>(position - 1);
        }
        else {
          start = 0;
        }
        end = static_cast<int>(std::min(position + alen + n_digits + 1, hlen));
        static constexpr int base_value {10};
        number = std::strtoll(std::next(header, position + alen), nullptr, base_value);

        return true;
      }

    return false;
  }


  auto find_abundance(struct seqinfo_s & seqinfo, struct Seq_stats & seq_stats, uint64_t lineno,
                      bool opt_usearch_abundance, int64_t opt_append_abundance) -> void
  {
    auto const & header_view = seqinfo.header_view;

    /* read size/abundance annotation */
    int64_t abundance = 0;
    int start = 0;
    int end = 0;
    int64_t number = 0;

    if (opt_usearch_abundance)
      {
        /* (^|;)size=([0-9]+)(;|$) */

        if (find_usearch_abundance(header_view, start, end, number))
          {
            if (number <= 0) {
              fatal(error_prefix, "Illegal abundance value on line ", lineno, ":\n",
                    header_view.data(), "\nAbundance values should be positive integers.");
            }
            abundance = number;
          }
      }
    else
      {
        /* (_)([0-9]+)$ */

        if (find_swarm_abundance(header_view, start, end, number))
          {
            if (number <= 0) {
              fatal(error_prefix, "Illegal abundance value on line ", lineno, ":\n",
                    header_view.data(), "\nAbundance values should be positive integers.");
            }
            abundance = number;
          }
      }

    if (abundance == 0)
      {
        start = static_cast<int>(header_view.size());
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
                seq_stats.missingabundance_header = header_view.data();
              }
          }
      }

    seqinfo.abundance = static_cast<uint64_t>(abundance);
    seqinfo.abundance_start = start;
    seqinfo.abundance_end = end;
  }


  auto abort_if_duplicated_identifier(struct seqinfo_s const * hdrfound,
                                      View<char const> const id_view) -> void {
    if (hdrfound == nullptr) { return; }
    std::string const id_str {id_view.data(), id_view.size()};
    fatal(error_prefix,
          "Duplicated sequence identifier: ",
          id_str);
  }


  auto abort_if_duplicated_sequences(struct Seq_stats const & seq_stats) -> void {
    if (not seq_stats.has_duplicates) { return; }
    fatal(error_prefix,
          "some fasta entries have identical sequences.\n"
          "Swarm expects dereplicated fasta files.\n"
          "Such files can be produced with swarm or vsearch:\n"
          " swarm -d 0 -w derep.fasta -o /dev/null input.fasta\n"
          "or\n"
          " vsearch --derep_fulllength input.fasta --sizein --sizeout --output derep.fasta");
  }


  auto abort_if_missing_abundance(struct Seq_stats const & seq_stats) -> void {
    if (seq_stats.missingabundance == 0) { return; }
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


  auto sort_index_if_need_be(struct Parameters const & parameters,
                             std::vector<struct seqinfo_s> & seqindex_v) -> void {
    struct Progress_status progress;
    progress_init(progress, "Abundance sorting:", 1, parameters);

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
      return lhs.header_view < rhs.header_view;
    };

    if (not std::is_sorted(seqindex_v.begin(), seqindex_v.end(),
                           compare_entries)) {
      std::sort(seqindex_v.begin(), seqindex_v.end(), compare_entries);
    }
    progress_done(progress);
  }


  auto print_user_report(struct Parameters const & parameters,
                         struct Seq_stats const & seq_stats) -> void {
    static_cast<void>(std::fprintf(parameters.logfile,
                                   "Database info:     %" PRIu64 " nt",
                                   seq_stats.nucleotides));
    static_cast<void>(std::fprintf(parameters.logfile,
                                   " in %u sequences,",
                                   seq_stats.n_sequences));
    static_cast<void>(std::fprintf(parameters.logfile,
                                   " longest %u nt\n",
                                   seq_stats.longest_sequence));
  }


  auto parse_fasta(struct Parameters const & parameters,
                   std::vector<char> & data_v) -> struct Parse_result
  {
    static constexpr unsigned int linealloc {2048};
    static constexpr unsigned int max_sequence_length {67108861};  // (2^26 - 3)
    // for longer sequences, 'zobrist_tab_byte_base' is bigger than 8 x
    // 2^32 (512 x max_sequence_length) and cannot be addressed with
    // uint32 pointers, which leads to a segmentation fault
    static constexpr unsigned int max_header_length {16777216 - 1};  // 2^24 minus 1

    auto const map_nt = make_nt_map();
    struct Parse_result result;
    auto & seq_stats = result.stats;
    auto & entries = result.entries;
    uint64_t datalen {0};

    /* open input file or stream */

    assert(parameters.input_filename.c_str() != nullptr);  // filename is set to '-' (stdin) by default

    auto const input_fp_handle = fopen_input(parameters.input_filename.c_str());
    if (not input_fp_handle)
      {
        fatal(error_prefix, "Unable to open input data file (", parameters.input_filename.c_str(), ").\n");
      }

    auto const file_info = get_file_info(input_fp_handle.get(), parameters);
    warn_if_file_is_not_regular(parameters, file_info.is_regular);

    /* allocate space */
    initial_allocation(data_v, file_info.filesize);

    uint64_t filepos = 0;

    Line_buffer line_buf{linealloc};

    auto lineno = 1U;


    struct Progress_status progress;
    progress_init(progress, "Reading sequences:", file_info.filesize, parameters);

    ssize_t linelen = xgetline(& line_buf.data, & line_buf.capacity, input_fp_handle.get());
    if (linelen < 0)
      {
        *line_buf.data = 0;
        linelen = 0;
      }
    filepos += static_cast<unsigned long int>(linelen);

    while (*line_buf.data != '\0')
      {
        /* read header */
        /* the header ends at a space, cr, lf or null character */

        if (*line_buf.data != '>') {
          fatal(error_prefix, "Illegal header line in fasta file.");
        }

        struct Entry entry;

        auto headerlen = static_cast<unsigned int>
          (std::strcspn(std::next(line_buf.data), " \r\n"));

        seq_stats.longestheader = std::max(headerlen, seq_stats.longestheader);

        if (seq_stats.longestheader > max_header_length) {
          fatal(error_prefix, "Headers longer than 16,777,215 symbols are not supported.");
        }

        /* store the line number */

        entry.lineno = lineno;


        /* store the header */

        linear_resize_if_need_be(data_v, datalen + headerlen + 1);
        std::memcpy(&data_v[datalen], std::next(line_buf.data), headerlen);
        data_v[datalen + headerlen] = '\0';
        entry.header.offset = datalen;
        entry.header.length = headerlen;  // '>' removed, so header is one byte shorter
        datalen += headerlen + 1;

        /* get next line */

        linelen = xgetline(& line_buf.data, & line_buf.capacity, input_fp_handle.get());
        if (linelen < 0)
          {
            *line_buf.data = '\0';
            linelen = 0;
          }
        filepos += static_cast<unsigned long int>(linelen);

        ++lineno;


        /* store a dummy sequence length */

        auto length = 0U;


        /* read and store sequence */

        uint64_t nt_buffer {0};
        auto nt_bufferlen = 0U;
        static constexpr unsigned int nt_buffersize {4 * sizeof(nt_buffer)};
        static constexpr unsigned char null_char = '\0';
        static constexpr int new_line {10};
        static constexpr int carriage_return {13};
        static constexpr int start_chars_range {32};  // visible ascii chars: 32-126
        static constexpr int end_chars_range {126};
        entry.sequence.offset = datalen;

        while ((*line_buf.data != 0) and (*line_buf.data != '>'))
          {
            auto character = null_char;
            auto * line_ptr = line_buf.data;
            while ((character = static_cast<unsigned char>(*line_ptr)) != null_char)
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
                        linear_resize_if_need_be(data_v, datalen + sizeof(nt_buffer));
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

            linelen = xgetline(& line_buf.data, & line_buf.capacity, input_fp_handle.get());
            if (linelen < 0)
              {
                *line_buf.data = 0;
                linelen = 0;
              }
            filepos += static_cast<unsigned long int>(linelen);

            ++lineno;
          }

        /* fill in real length */

        entry.sequence.length = length;

        if (length == 0)
          {
            fatal(error_prefix, "Empty sequence found on line ", lineno - 1, ".");
          }

        seq_stats.nucleotides += length;
        seq_stats.longest_sequence = std::max(length, seq_stats.longest_sequence);


        /* save remaining padded 64-bit value with nt's, if any */

        if (nt_bufferlen > 0)
          {
            linear_resize_if_need_be(data_v, datalen + sizeof(nt_buffer));
            std::memcpy(&data_v[datalen], & nt_buffer, sizeof(nt_buffer));
            datalen += sizeof(nt_buffer);

            nt_buffer = 0;
            nt_bufferlen = 0;  // that value is never read again, all tests pass without it
          }

        ++seq_stats.n_sequences;
        entries.push_back(entry);

        if (file_info.is_regular) {
          progress_update(progress, filepos);
        }
      }
    progress_done(progress);

    // Line_buffer is destroyed on return; indexing/hashing in
    // build_index can use the released memory
    return result;
  }


  auto build_index(struct Parameters const & parameters,
                   Zobrist const & zobrist,
                   std::vector<char> & data_v,
                   std::vector<struct Entry> const & entries,
                   struct Seq_stats & seq_stats,
                   std::vector<struct seqinfo_s> & seqindex_v) -> void
  {
    /* set up hash to check for unique headers */

    const uint64_t hdrhashsize {2ULL * seq_stats.n_sequences};
    std::vector<struct seqinfo_s *> hdrhashtable(hdrhashsize);

    /* set up hash to check for unique sequences */

    const uint64_t seqhashsize {2ULL * seq_stats.n_sequences};

    std::vector<struct seqinfo_s *> seqhashtable;

    if (parameters.opt_differences > 1) {
      seqhashtable.resize(seqhashsize);
    }

    /* create indices */

    seqindex_v.resize(seq_stats.n_sequences);

    struct Progress_status progress_idx;
    progress_init(progress_idx, "Indexing database:", seq_stats.n_sequences, parameters);
    auto counter = 0ULL;
    for (auto & a_sequence: seqindex_v) {

        /* get header */
        a_sequence.header_view = View<char const>{
          &data_v[entries[counter].header.offset],
          entries[counter].header.length};

        /* and sequence */
        const auto seqlen = static_cast<unsigned int>(entries[counter].sequence.length);
        a_sequence.seqlen = seqlen;
        a_sequence.seq = &data_v[entries[counter].sequence.offset];

        /* get amplicon abundance */
        find_abundance(a_sequence, seq_stats, entries[counter].lineno, parameters.opt_usearch_abundance, parameters.opt_append_abundance);

        auto const headerlen_signed = static_cast<int>(a_sequence.header_view.size());
        if ((a_sequence.abundance_start == 0) and
            (a_sequence.abundance_end == headerlen_signed)) {
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
            id_len = headerlen_signed - a_sequence.abundance_end;
          }

        auto const id_view = a_sequence.header_view.subview(
          static_cast<std::size_t>(id_start),
          static_cast<std::size_t>(id_len));

        const auto hdrhash = zobrist.hash(id_view.data(),
                                          4 * static_cast<unsigned int>(id_len));

        a_sequence.hdrhash = hdrhash;
        uint64_t hdrhashindex = hdrhash % hdrhashsize;

        struct seqinfo_s const * hdrfound {nullptr};

        while ((hdrfound = hdrhashtable[hdrhashindex]) != nullptr)
          {
            if (hdrfound->hdrhash == hdrhash)
              {
                int hit_id_start {0};
                int hit_id_len {0};

                auto const hit_headerlen_signed =
                  static_cast<int>(hdrfound->header_view.size());
                if (hdrfound->abundance_start > 0)
                  {
                    hit_id_start = 0;
                    hit_id_len = hdrfound->abundance_start;
                  }
                else
                  {
                    hit_id_start = hdrfound->abundance_end;
                    hit_id_len = hit_headerlen_signed - hdrfound->abundance_end;
                  }

                auto const hit_id_view = hdrfound->header_view.subview(
                  static_cast<std::size_t>(hit_id_start),
                  static_cast<std::size_t>(hit_id_len));
                if (id_view == hit_id_view) {
                  break;
                }
              }

            hdrhashindex = (hdrhashindex + 1) % hdrhashsize;
          }

        abort_if_duplicated_identifier(hdrfound, id_view);

        hdrhashtable[hdrhashindex] = &a_sequence;

        /* hash sequence */
        a_sequence.seqhash = zobrist.hash(a_sequence.seq, a_sequence.seqlen);

        if (parameters.opt_differences > 1)
          {
            // refactoring: extract to a free function (not trivial)
            /* Check for duplicated sequences using hash table,  */
            /* but only for d > 1. Handled internally for d = 1. */

            uint64_t seqhashindex = a_sequence.seqhash % seqhashsize;
            struct seqinfo_s const * seqfound {nullptr};

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
                seq_stats.has_duplicates = true;
                break;
              }
            seqhashtable[seqhashindex] = &a_sequence;
          }

        progress_update(progress_idx, counter);
        ++counter;
      }

    abort_if_duplicated_sequences(seq_stats);

    progress_done(progress_idx);

    abort_if_missing_abundance(seq_stats);
    sort_index_if_need_be(parameters, seqindex_v);
    print_user_report(parameters, seq_stats);
  }

} // end of anonymous namespace


// ----- class Data -----

Data::Data(struct Parameters const & parameters)
{
  auto parse_result = parse_fasta(parameters, data_);

  // Construct the Zobrist tables now that parse_fasta has determined
  // the longest header and sequence. The +2 budgets two insertions
  // for the variant enumeration in variants.cc.
  auto const & stats = parse_result.stats;
  longest_ = stats.longest_sequence;
  auto const zobrist_len = std::max(4 * stats.longestheader, stats.longest_sequence + 2);
  zobrist_p_.reset(new Zobrist(zobrist_len));

  build_index(parameters, *zobrist_p_, data_, parse_result.entries, parse_result.stats, seqindex_);
}


auto Data::info(uint64_t const seqno) const -> struct seqinfo_s const &
{
  assert(not seqindex_.empty());  // db_read() / Data ctor must run first
  // bound-check is redundant with -D_GLIBCXX_DEBUG (operator[] is
  // already checked under libstdc++ debug mode), kept here so the
  // precondition is enforced in any assert-enabled build
  assert(seqno < seqindex_.size());
  return seqindex_[seqno];
}


auto Data::sequence(uint64_t const seqno) const -> char const *
{
  return info(seqno).seq;
}


auto Data::sequence_length(uint64_t const seqno) const -> unsigned int
{
  return info(seqno).seqlen;
}


auto Data::sequence_hash(uint64_t const seqno) const -> uint64_t
{
  return info(seqno).seqhash;
}


auto Data::header_view(uint64_t const seqno) const -> View<char const>
{
  return info(seqno).header_view;
}


auto Data::abundance(uint64_t const seqno) const -> uint64_t
{
  return info(seqno).abundance;
}


// refactoring: decompress sequence (4 nt at a time)
// - need a const vector<string> byte_decode = { "AAAA", "AAAC", "AAAG", ...
// - need a std::string buffer of capacity = length + 3 + 1,
// - for each compressed byte in the view of length (len + 3) % 4 bytes,
//   for (auto const compressed_byte : compressed_bytes) {
//       auto const s_compressed_byte = static_cast<unsigned char>(compressed_byte);
//       buffer += byte_decode[s_compressed_byte];
//   buffer[len] = '\0';  //
//   std::fprintf(fastaout_fp, "%.*s\n", len, buffer.c_str());
// benchmarck to check which way is faster
auto Data::fprintseq(std::FILE * stream, unsigned int const seqno) const -> void
{
  static constexpr std::array<char, 32> sym_nt =
    {'-', 'A', 'C', 'G', 'T', ' ', ' ', ' ',
     ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
     ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ',
     ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '};
  auto const len = sequence_length(seqno);
  auto const * const seqptr = sequence(seqno);
  static std::vector<char> buffer(longest_sequence() + 1, '\0');

  // decode to nucleotides (A, C, G and T)
  for (auto i = 0U; i < len; ++i) {
    buffer[i] = sym_nt[1 + nt_extract(seqptr, i)];
  }
  buffer[len] = '\0';

  std::fprintf(stream, "%.*s\n", len, buffer.data());
}


auto Data::fprint_id(std::FILE * stream, uint64_t const seqno,
                     bool const opt_usearch_abundance,
                     int64_t const opt_append_abundance) const -> void
{
  auto const & seqinfo = info(seqno);
  auto const * hdrstr = seqinfo.header_view.data();
  auto const hdrlen = static_cast<int>(seqinfo.header_view.size());
  auto const abundance_value = seqinfo.abundance;

  // if abundance is missing and if user says that a missing abundance is ok, then...
  if ((opt_append_abundance != 0) and (seqinfo.abundance_start == seqinfo.abundance_end)) {
    if (opt_usearch_abundance) {
      std::fprintf(stream, "%.*s;size=%" PRIu64 ";", hdrlen, hdrstr, abundance_value);
    }
    else {
      std::fprintf(stream, "%.*s_%" PRIu64, hdrlen, hdrstr, abundance_value);
    }
  }
  else {
    std::fprintf(stream, "%.*s", hdrlen, hdrstr);
  }
}


auto Data::fprint_id_noabundance(std::FILE * stream, uint64_t const seqno,
                                 bool const opt_usearch_abundance) const -> void
{
  auto const & seqinfo = info(seqno);
  auto const * hdrstr = seqinfo.header_view.data();
  auto const hdrlen = static_cast<int>(seqinfo.header_view.size());
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


auto Data::fprint_id_with_new_abundance(std::FILE * stream,
                                        uint64_t const seqno,
                                        uint64_t const new_abundance,
                                        bool const opt_usearch_abundance) const -> void
{
  auto const & seqinfo = info(seqno);

  auto const * const hdrstr = seqinfo.header_view.data();
  auto const hdrlen = static_cast<int>(seqinfo.header_view.size());

  if (opt_usearch_abundance) {
    std::fprintf(stream,
                 "%.*s%ssize=%" PRIu64 ";%.*s",
                 seqinfo.abundance_start,
                 hdrstr,
                 seqinfo.abundance_start > 0 ? ";" : "",
                 new_abundance,
                 hdrlen - seqinfo.abundance_end,
                 std::next(hdrstr, seqinfo.abundance_end));
  }
  else {
    std::fprintf(stream,
                 "%.*s_%" PRIu64,
                 seqinfo.abundance_start,
                 hdrstr,
                 new_abundance);
  }
}


