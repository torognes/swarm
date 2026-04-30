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
#include "utils/hasher_fnv1a.h"
#include "utils/hasher_generic.h"
#include "utils/input_output.h"
#include "utils/nt_codec.h"
#include "utils/progress.h"
#include "utils/seq_index.h"
#include "utils/view.h"
#include "utils/xgetline.h"
#include <algorithm>  // std::all_of() std::copy_n() std::find() std::find_if_not() std::max() std::min() std::search() std::sort()
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
#include <unordered_set>
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
  constexpr unsigned int max_header_length {16777216 - 1};  // 2^24 minus 1
  constexpr unsigned int max_sequence_length {67108861};  // (2^26 - 3)
  // for longer sequences, 'zobrist_tab_byte_base' is bigger than 8 x
  // 2^32 (512 x max_sequence_length) and cannot be addressed with
  // uint32 pointers, which leads to a segmentation fault

  // Nucleotide character classification: the lookup table built by
  // make_nt_classifier() returns one of these for every ASCII byte.
  // The four nucleotide values are also the packed 2-bit encoding,
  // so they can be passed straight to Nt_packer::push(). Ordering
  // matters: nucleotides are < nt_class_skip < nt_class_illegal so
  // the hot-path test is a single comparison.
  constexpr uint8_t nt_class_a       {0};
  constexpr uint8_t nt_class_c       {1};
  constexpr uint8_t nt_class_g       {2};
  constexpr uint8_t nt_class_t       {3};
  constexpr uint8_t nt_class_skip    {4};
  constexpr uint8_t nt_class_illegal {5};

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


  // Result of a successful abundance-annotation parse. 'found' is the
  // discriminator: when false, the other fields are meaningless.
  struct Abundance_match {
    int     start  {0};
    int     end    {0};
    int64_t number {0};
    bool    found  {false};
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


  auto make_nt_classifier() -> std::array<uint8_t, n_chars> {
    // every ascii byte falls into exactly one of: nucleotide (A/C/G/T/U,
    // case insensitive) -> packed 2-bit encoding; line terminator
    // (CR or LF) -> silently skipped; anything else -> fatal error
    std::array<uint8_t, n_chars> table;
    table.fill(nt_class_illegal);
    table['A'] = nt_class_a;  table['a'] = nt_class_a;
    table['C'] = nt_class_c;  table['c'] = nt_class_c;
    table['G'] = nt_class_g;  table['g'] = nt_class_g;
    table['T'] = nt_class_t;  table['t'] = nt_class_t;
    table['U'] = nt_class_t;  table['u'] = nt_class_t;
    table['\n'] = nt_class_skip;
    table['\r'] = nt_class_skip;
    return table;
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


  // Read one line into line_buf and bump filepos by the number of bytes
  // consumed. On read failure, leave the buffer empty (first byte set
  // to '\0') so callers can use the same end-of-input sentinel.
  auto read_next_line(Line_buffer & line_buf, std::FILE * stream,
                      uint64_t & filepos) -> void
  {
    auto const linelen = xgetline(& line_buf.data, & line_buf.capacity, stream);
    if (linelen < 0) {
      *line_buf.data = '\0';
      return;
    }
    filepos += static_cast<unsigned long int>(linelen);
  }


  // Pack 4 nucleotides per byte into a 64-bit accumulator and flush
  // it to data_v as a fixed-size memcpy whenever it fills up. A final
  // flush() at end-of-sequence writes the partially-filled buffer
  // padded with zeros (so the on-disk layout is unchanged).
  struct Nt_packer {
    uint64_t buffer {0};
    unsigned int filled {0};
    static constexpr unsigned int capacity {4 * sizeof(buffer)};  // 32 bases per uint64

    auto push(uint64_t const mapped_minus_one,
              std::vector<char> & data_v, uint64_t & datalen) -> void
    {
      buffer |= mapped_minus_one << (2 * filled);
      ++filled;
      if (filled == capacity) { flush(data_v, datalen); }
    }

    auto flush(std::vector<char> & data_v, uint64_t & datalen) -> void
    {
      linear_resize_if_need_be(data_v, datalen + sizeof(buffer));
      std::memcpy(&data_v[datalen], &buffer, sizeof(buffer));
      datalen += sizeof(buffer);
      buffer = 0;
      filled = 0;
    }
  };


  // Validate the '>' header line, copy the header bytes (everything
  // after '>' up to the first space, CR or LF) into data_v, and fill
  // entry.header. Updates seq_stats.longestheader and aborts when
  // max_header_length is exceeded.
  auto store_header(Line_buffer const & line_buf,
                    struct Entry & entry,
                    std::vector<char> & data_v,
                    uint64_t & datalen,
                    struct Seq_stats & seq_stats) -> void
  {
    if (*line_buf.data != '>') {
      fatal(error_prefix, "Illegal header line in fasta file.");
    }

    auto const headerlen = static_cast<unsigned int>
      (std::strcspn(std::next(line_buf.data), " \r\n"));

    seq_stats.longestheader = std::max(headerlen, seq_stats.longestheader);

    if (seq_stats.longestheader > max_header_length) {
      fatal(error_prefix, "Headers longer than 16,777,215 symbols are not supported.");
    }

    linear_resize_if_need_be(data_v, datalen + headerlen + 1);
    std::copy_n(std::next(line_buf.data), headerlen, &data_v[datalen]);
    data_v[datalen + headerlen] = '\0';
    entry.header.offset = datalen;
    entry.header.length = headerlen;  // '>' removed, so header is one byte shorter
    datalen += headerlen + 1;
  }


  // Read sequence lines starting at line_buf (already loaded with the
  // first line after a header) until the next header or end of input.
  // Pack nucleotides into data_v via Nt_packer, validate characters,
  // enforce max_sequence_length, and write entry.sequence + the
  // running counters in seq_stats. Stops with line_buf holding the
  // line that broke the loop ('>' or '\0').
  auto parse_sequence_body(Line_buffer & line_buf, std::FILE * stream,
                           std::array<uint8_t, n_chars> const & classify,
                           std::vector<char> & data_v, uint64_t & datalen,
                           uint64_t & filepos, unsigned int & lineno,
                           struct Entry & entry,
                           struct Seq_stats & seq_stats) -> void
  {
    static constexpr unsigned char null_char = '\0';
    static constexpr int start_chars_range {32};  // visible ascii chars: 32-126
    static constexpr int end_chars_range {126};

    Nt_packer packer;
    auto length = 0U;
    entry.sequence.offset = datalen;

    while ((*line_buf.data != 0) and (*line_buf.data != '>'))
      {
        auto * line_ptr = line_buf.data;
        unsigned char character {};
        while ((character = static_cast<unsigned char>(*line_ptr)) != null_char)
          {
            line_ptr = std::next(line_ptr);
            auto const category = classify[character];
            if (category < nt_class_skip)
              {
                packer.push(category, data_v, datalen);
                ++length;
              }
            else if (category == nt_class_illegal)
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
            // else: nt_class_skip (CR or LF), silently ignored
          }

        /* check length of longest sequence */
        if (length > max_sequence_length) {
          fatal(error_prefix, "Sequences longer than 67,108,861 symbols are not supported.");
        }

        read_next_line(line_buf, stream, filepos);

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

    if (packer.filled > 0) {
      packer.flush(data_v, datalen);
    }
  }


  auto find_swarm_abundance(View<char> const header_view) -> Abundance_match
  {
    /*
      Identify the first occurence of the pattern (_)([0-9]+)$
      in the header string.
    */

    static constexpr std::size_t max_digits {20};  // 20 digits at most (abundance > 10^20)

    auto const is_digit = [](char const character) noexcept -> bool {
      return (character >= '0') and (character <= '9');
    };

    // Find the last '_' via reverse scan over the header view.
    auto const r_underscore = std::find(header_view.crbegin(),
                                        header_view.crend(), '_');
    if (r_underscore == header_view.crend()) {
      return Abundance_match{};
    }

    // base() of a reverse iterator points one past the matched element,
    // i.e. at the first byte after the '_'.
    auto const * const digits_begin = r_underscore.base();
    auto const * const digits_end   = header_view.cend();
    auto const n_digits = static_cast<std::size_t>(
      std::distance(digits_begin, digits_end));

    if ((n_digits == 0) or (n_digits > max_digits)) {
      return Abundance_match{};
    }
    if (not std::all_of(digits_begin, digits_end, is_digit)) {
      return Abundance_match{};
    }

    auto const underscore_offset = std::distance(header_view.cbegin(),
                                                 std::prev(digits_begin));
    assert(underscore_offset >= 0);
    assert(underscore_offset <= std::numeric_limits<int>::max());
    assert(n_digits <= static_cast<std::size_t>(std::numeric_limits<int>::max()));

    Abundance_match match;
    match.start = static_cast<int>(underscore_offset);
    match.end   = match.start + 1 + static_cast<int>(n_digits);

    // strtoll still requires null-termination at the end of the digit run;
    // header_view points into Data::data_, where each header is followed
    // by a '\0' byte written at parse time.
    // refactoring: capture strtoll's end pointer and check errno == ERANGE
    // to detect overflow (n_digits is bounded above by max_digits = 20,
    // which can exceed int64_t's 19-digit range).
    static constexpr int base_value {10};
    match.number = std::strtoll(digits_begin, nullptr, base_value);
    match.found  = true;
    return match;
  }


  auto find_usearch_abundance(View<char> const header_view) -> Abundance_match
  {
    /*
      Identify the first occurence of the pattern (^|;)size=([0-9]+)(;|$)
      in the header string.
    */

    static constexpr char attribute[] {"size="};
    static constexpr std::size_t alen {sizeof(attribute) - 1};  // exclude trailing '\0'

    auto const is_digit = [](char const character) noexcept -> bool {
      return (character >= '0') and (character <= '9');
    };

    auto const * const header_begin = header_view.cbegin();
    auto const * const header_end   = header_view.cend();
    auto const * search_from = header_begin;

    while (search_from != header_end)
      {
        auto const * const match = std::search(search_from, header_end,
                                               std::begin(attribute),
                                               std::next(std::begin(attribute), alen));
        if (match == header_end) {
          return Abundance_match{};
        }

        auto const * const digits_begin = std::next(match, alen);

        /* left context: start of header or ';' */
        bool const left_ok = (match == header_begin)
                          or (*std::prev(match) == ';');

        /* digit run, then right context: end of header or ';' */
        auto const * const digits_end = std::find_if_not(digits_begin, header_end,
                                                         is_digit);
        auto const n_digits = std::distance(digits_begin, digits_end);
        bool const right_ok = (digits_end == header_end)
                           or (*digits_end == ';');

        if (left_ok and (n_digits > 0) and right_ok)
          {
            auto const match_offset = std::distance(header_begin, match);
            assert(match_offset >= 0);
            assert(match_offset <= std::numeric_limits<int>::max());

            Abundance_match result;
            result.start = (match_offset > 0) ? static_cast<int>(match_offset - 1) : 0;

            // include the trailing ';' when present, otherwise stop at end
            auto end_offset = std::distance(header_begin, digits_end);
            if (digits_end != header_end) {
              ++end_offset;
            }
            assert(end_offset <= std::numeric_limits<int>::max());
            result.end = static_cast<int>(end_offset);

            // strtoll still requires null-termination at the end of the
            // digit run; the digit run is always followed by either ';'
            // or the '\0' at the end of the header in Data::data_.
            static constexpr int base_value {10};
            result.number = std::strtoll(digits_begin, nullptr, base_value);
            result.found  = true;
            return result;
          }

        // skip past this 'size=' and keep scanning
        search_from = digits_begin;
      }

    return Abundance_match{};
  }


  auto find_abundance(struct seqinfo_s & seqinfo, struct Seq_stats & seq_stats, uint64_t lineno,
                      bool opt_usearch_abundance, int64_t opt_append_abundance) -> void
  {
    auto const & header_view = seqinfo.header_view;

    /* read size/abundance annotation */
    auto const match = opt_usearch_abundance
      ? find_usearch_abundance(header_view)  /* (^|;)size=([0-9]+)(;|$) */
      : find_swarm_abundance(header_view);   /* (_)([0-9]+)$ */

    int64_t abundance = 0;
    int start = match.start;
    int end   = match.end;

    if (match.found)
      {
        if (match.number <= 0) {
          fatal(error_prefix, "Illegal abundance value on line ", lineno, ":\n",
                header_view.data(), "\nAbundance values should be positive integers.");
        }
        abundance = match.number;
      }
    else
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
    Progress progress("Abundance sorting:", 1, parameters);

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
    progress.done();
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

    auto const classify = make_nt_classifier();
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


    Progress progress("Reading sequences:", file_info.filesize, parameters);

    read_next_line(line_buf, input_fp_handle.get(), filepos);

    while (*line_buf.data != '\0')
      {
        /* read header */
        /* the header ends at a space, cr, lf or null character */

        struct Entry entry;
        entry.lineno = lineno;

        store_header(line_buf, entry, data_v, datalen, seq_stats);

        /* get next line */

        read_next_line(line_buf, input_fp_handle.get(), filepos);

        ++lineno;


        /* read and store sequence */

        parse_sequence_body(line_buf, input_fp_handle.get(), classify,
                            data_v, datalen, filepos, lineno,
                            entry, seq_stats);

        ++seq_stats.n_sequences;
        entries.push_back(entry);

        if (file_info.is_regular) {
          progress.update(filepos);
        }
      }
    progress.done();

    // Line_buffer is destroyed on return; indexing/hashing in
    // build_index can use the released memory
    return result;
  }


  // Populate header_view and the (seq, seqlen) pointer pair from a
  // parsed Entry into its destination seqinfo slot.
  auto populate_views_from_entry(struct seqinfo_s & a_sequence,
                                 struct Entry const & entry,
                                 std::vector<char> & data_v) -> void
  {
    a_sequence.header_view = View<char>{
      &data_v[entry.header.offset],
      entry.header.length};
    a_sequence.seqlen = static_cast<unsigned int>(entry.sequence.length);
    a_sequence.seq    = &data_v[entry.sequence.offset];
  }


  // Compute the identifier subview within header_view, given the
  // abundance range already filled in by find_abundance(). Aborts if
  // the identifier would be empty.
  auto compute_identifier_view(struct seqinfo_s const & a_sequence) -> View<char>
  {
    auto const headerlen_signed = static_cast<int>(a_sequence.header_view.size());
    if ((a_sequence.abundance_start == 0) and
        (a_sequence.abundance_end == headerlen_signed)) {
      fatal(error_prefix, "Empty sequence identifier.");
    }

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

    return a_sequence.header_view.subview(
      static_cast<std::size_t>(id_start),
      static_cast<std::size_t>(id_len));
  }


  // Insert id_view into the dedup set; abort with the offending
  // identifier in the error message if it was already present.
  auto register_unique_identifier(
    std::unordered_set<View<char>, GenericHash<fnv1a>> & seen_identifiers,
    View<char> const id_view) -> void
  {
    auto const insertion = seen_identifiers.insert(id_view);
    if (not insertion.second) {
      std::string const id_str {id_view.data(), id_view.size()};
      fatal(error_prefix, "Duplicated sequence identifier: ", id_str);
    }
  }


  // Open-addressed lookup over seqhashtable, used when d > 0
  // (d = 0 is the dereplication mode and accepts duplicates).
  // Returns true if an identical sequence was already inserted;
  // otherwise records a_sequence at the probed slot and returns false.
  auto is_duplicate_sequence(std::vector<struct seqinfo_s *> & seqhashtable,
                             uint64_t const seqhashsize,
                             struct seqinfo_s & a_sequence) -> bool
  {
    uint64_t seqhashindex = a_sequence.seqhash % seqhashsize;
    struct seqinfo_s const * seqfound {nullptr};

    while ((seqfound = seqhashtable[seqhashindex]) != nullptr)
      {
        if ((seqfound->seqhash == a_sequence.seqhash) and
            (seqfound->seqlen == a_sequence.seqlen) and
            std::equal(seqfound->seq,
                       std::next(seqfound->seq, nt_bytelength(a_sequence.seqlen)),
                       a_sequence.seq)) {
          return true;
        }
        seqhashindex = (seqhashindex + 1) % seqhashsize;
      }

    seqhashtable[seqhashindex] = &a_sequence;
    return false;
  }


  // Pass 1: header-side work for every entry. Populates views,
  // extracts the abundance annotation, and aborts on duplicated or
  // empty identifiers. The seen_identifiers set is local so its
  // buckets are released before the sequence pass allocates its own
  // hashtable.
  auto index_headers(struct Parameters const & parameters,
                     std::vector<char> & data_v,
                     std::vector<struct Entry> const & entries,
                     struct Seq_stats & seq_stats,
                     std::vector<struct seqinfo_s> & seqindex_v,
                     Progress & progress_idx) -> void
  {
    std::unordered_set<View<char>, GenericHash<fnv1a>> seen_identifiers;
    seen_identifiers.reserve(seq_stats.n_sequences);

    auto counter = 0ULL;
    for (auto & a_sequence: seqindex_v) {
        populate_views_from_entry(a_sequence, entries[counter], data_v);

        /* get amplicon abundance */
        find_abundance(a_sequence, seq_stats, entries[counter].lineno,
                       parameters.opt_usearch_abundance, parameters.opt_append_abundance);

        auto const id_view = compute_identifier_view(a_sequence);
        register_unique_identifier(seen_identifiers, id_view);

        progress_idx.update(counter);
        ++counter;
      }
  }


  // Pass 2: sequence-side work for every entry. Computes the Zobrist
  // hash and, when d > 0, checks for duplicated sequences (d = 0 is
  // the dereplication mode and accepts duplicates). Stops at the first
  // duplicate; the caller calls abort_if_duplicated_sequences() afterwards.
  auto index_sequences(struct Parameters const & parameters,
                       Zobrist const & zobrist,
                       struct Seq_stats & seq_stats,
                       std::vector<struct seqinfo_s> & seqindex_v,
                       Progress & progress_idx) -> void
  {
    const uint64_t seqhashsize {2ULL * seq_stats.n_sequences};
    std::vector<struct seqinfo_s *> seqhashtable;
    if (parameters.opt_differences > 0) {
      seqhashtable.resize(seqhashsize);
    }

    auto counter = 0ULL;
    for (auto & a_sequence: seqindex_v) {
        a_sequence.seqhash = zobrist.hash(a_sequence.seq, a_sequence.seqlen);

        if ((parameters.opt_differences > 0) and
            is_duplicate_sequence(seqhashtable, seqhashsize, a_sequence))
          {
            seq_stats.has_duplicates = true;
            break;
          }

        progress_idx.update(seq_stats.n_sequences + counter);
        ++counter;
      }
  }


  auto build_index(struct Parameters const & parameters,
                   Zobrist const & zobrist,
                   std::vector<char> & data_v,
                   std::vector<struct Entry> const & entries,
                   struct Seq_stats & seq_stats,
                   std::vector<struct seqinfo_s> & seqindex_v) -> void
  {
    seqindex_v.resize(seq_stats.n_sequences);

    // One progress bar drives both passes: pass 1 contributes the
    // first half of the count, pass 2 the second half.
    Progress progress_idx("Indexing database:",
                          2ULL * seq_stats.n_sequences, parameters);

    index_headers(parameters, data_v, entries, seq_stats, seqindex_v, progress_idx);
    index_sequences(parameters, zobrist, seq_stats, seqindex_v, progress_idx);

    abort_if_duplicated_sequences(seq_stats);

    progress_idx.done();

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


auto Data::sequence_view(uint64_t const seqno) const -> Sequence
{
  auto const & rec = info(seqno);
  return {View<char>{rec.seq, nt_bytelength(rec.seqlen)}, rec.seqlen};
}


auto Data::sequence_hash(uint64_t const seqno) const -> uint64_t
{
  return info(seqno).seqhash;
}


auto Data::header_view(uint64_t const seqno) const -> View<char>
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
  auto const seq = sequence_view(seqno);
  static std::vector<char> buffer(longest_sequence() + 1, '\0');

  // decode to nucleotides (A, C, G and T)
  for (auto i = 0U; i < seq.length; ++i) {
    buffer[i] = sym_nt[1 + nt_extract(seq.encoded.data(), i)];
  }
  buffer[seq.length] = '\0';

  std::fprintf(stream, "%.*s\n", seq.length, buffer.data());
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


