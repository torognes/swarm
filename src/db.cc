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
#include "utils/line_buffer.h"
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
#include <vector>


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
  // so they can be passed straight to Nt_packer::push() after a cast
  // to the underlying type. Ordering matters: bases < skip < illegal,
  // so the hot-path test is a single "<" comparison.
  enum struct Nt_class : uint8_t {
    a       = 0,
    c       = 1,
    g       = 2,
    t       = 3,
    skip    = 4,
    illegal = 5,
  };

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


  auto make_nt_classifier() -> std::array<Nt_class, n_chars> {
    // every ascii byte falls into exactly one of: nucleotide (A/C/G/T/U,
    // case insensitive) -> packed 2-bit encoding; line terminator
    // (CR or LF) -> silently skipped; anything else -> fatal error
    std::array<Nt_class, n_chars> table;
    table.fill(Nt_class::illegal);
    table['A'] = Nt_class::a;  table['a'] = Nt_class::a;
    table['C'] = Nt_class::c;  table['c'] = Nt_class::c;
    table['G'] = Nt_class::g;  table['g'] = Nt_class::g;
    table['T'] = Nt_class::t;  table['t'] = Nt_class::t;
    table['U'] = Nt_class::t;  table['u'] = Nt_class::t;
    table['\n'] = Nt_class::skip;
    table['\r'] = Nt_class::skip;
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
      fatal("Unable to fstat on input file (", parameters.input_filename.c_str(), ").\n");
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
    if (line_buf.peek_first() != '>') {
      fatal("Illegal header line in fasta file.");
    }

    auto const headerlen = static_cast<unsigned int>
      (std::strcspn(std::next(line_buf.data()), " \r\n"));

    seq_stats.longestheader = std::max(headerlen, seq_stats.longestheader);

    if (seq_stats.longestheader > max_header_length) {
      fatal("Headers longer than 16,777,215 symbols are not supported.");
    }

    linear_resize_if_need_be(data_v, datalen + headerlen + 1);
    std::copy_n(std::next(line_buf.data()), headerlen, &data_v[datalen]);
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
                           std::array<Nt_class, n_chars> const & classify,
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

    while ((not line_buf.empty()) and (line_buf.peek_first() != '>')) {
        auto const * line_ptr = line_buf.data();
        unsigned char character {};
        while ((character = static_cast<unsigned char>(*line_ptr)) != null_char) {
            line_ptr = std::next(line_ptr);
            auto const category = classify[character];
            if (category < Nt_class::skip) {
                packer.push(static_cast<uint8_t>(category), data_v, datalen);
                ++length;
              }
            else if (category == Nt_class::illegal) {
                if ((character >= start_chars_range) and (character <= end_chars_range)) {
                  fatal("Illegal character '", static_cast<char>(character),
                        "' in sequence on line ", lineno, ".");
                }
                else {
                  fatal("Illegal character (ascii no ", static_cast<unsigned int>(character),
                        ") in sequence on line ", lineno, ".");
                }
              }
            // else: Nt_class::skip (CR or LF), silently ignored
          }

        /* check length of longest sequence */
        if (length > max_sequence_length) {
          fatal("Sequences longer than 67,108,861 symbols are not supported.");
        }

        line_buf.read_next(stream, filepos);

        ++lineno;
      }

    /* fill in real length */

    entry.sequence.length = length;

    if (length == 0) {
        fatal("Empty sequence found on line ", lineno - 1, ".");
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

    static constexpr std::array<char, 5> attribute {{'s', 'i', 'z', 'e', '='}};

    auto const is_digit = [](char const character) noexcept -> bool {
      return (character >= '0') and (character <= '9');
    };

    auto const * const header_begin = header_view.cbegin();
    auto const * const header_end   = header_view.cend();
    auto const * search_from = header_begin;

    while (search_from != header_end) {
        auto const * const match = std::search(search_from, header_end,
                                               attribute.cbegin(),
                                               attribute.cend());
        if (match == header_end) {
          return Abundance_match{};
        }

        auto const * const digits_begin = std::next(match, attribute.size());

        /* left context: start of header or ';' */
        bool const left_ok = (match == header_begin)
                          or (*std::prev(match) == ';');

        /* digit run, then right context: end of header or ';' */
        auto const * const digits_end = std::find_if_not(digits_begin, header_end,
                                                         is_digit);
        auto const n_digits = std::distance(digits_begin, digits_end);
        bool const right_ok = (digits_end == header_end)
                           or (*digits_end == ';');

        if (left_ok and (n_digits > 0) and right_ok) {
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

    if (match.found) {
        if (match.number <= 0) {
          fatal("Illegal abundance value on line ", lineno, ":\n",
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
            if (seq_stats.missingabundance == 1) {
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
    fatal(
          "some fasta entries have identical sequences.\n"
          "Swarm expects dereplicated fasta files.\n"
          "Such files can be produced with swarm or vsearch:\n"
          " swarm -d 0 -w derep.fasta -o /dev/null input.fasta\n"
          "or\n"
          " vsearch --derep_fulllength input.fasta --sizein --sizeout --output derep.fasta");
  }


  auto abort_if_missing_abundance(struct Seq_stats const & seq_stats) -> void {
    if (seq_stats.missingabundance == 0) { return; }
    fatal("Abundance annotations not found for ",
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
    Progress const progress("Abundance sorting:", 1, parameters);

    auto compare_entries = [](struct seqinfo_s const & lhs,
                              struct seqinfo_s const & rhs) -> bool {
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
                   std::vector<char> & data_v) -> struct Parse_result {
    static constexpr unsigned int linealloc {2048};

    auto const classify = make_nt_classifier();
    struct Parse_result result;
    auto & seq_stats = result.stats;
    auto & entries = result.entries;
    uint64_t datalen {0};

    /* open input file or stream */

    assert(parameters.input_filename.c_str() != nullptr);  // filename is set to '-' (stdin) by default

    auto const input_fp_handle = fopen_input(parameters.input_filename.c_str());
    if (not input_fp_handle) {
        fatal("Unable to open input data file (", parameters.input_filename.c_str(), ").\n");
      }

    auto const file_info = get_file_info(input_fp_handle.get(), parameters);
    warn_if_file_is_not_regular(parameters, file_info.is_regular);

    /* allocate space */
    initial_allocation(data_v, file_info.filesize);

    uint64_t filepos = 0;

    Line_buffer line_buf{linealloc};

    auto lineno = 1U;


    Progress progress("Reading sequences:", file_info.filesize, parameters);

    line_buf.read_next(input_fp_handle.get(), filepos);

    while (not line_buf.empty()) {
        /* read header */
        /* the header ends at a space, cr, lf or null character */

        struct Entry entry;
        entry.lineno = lineno;

        store_header(line_buf, entry, data_v, datalen, seq_stats);

        /* get next line */

        line_buf.read_next(input_fp_handle.get(), filepos);

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
                                 std::vector<char> const & data_v) -> void {
    a_sequence.header_view = View<char>{
      &data_v[entry.header.offset],
      entry.header.length};
    a_sequence.seqlen = static_cast<unsigned int>(entry.sequence.length);
    a_sequence.seq    = &data_v[entry.sequence.offset];
  }


  // Compute the identifier subview within header_view, given the
  // abundance range already filled in by find_abundance(). Aborts if
  // the identifier would be empty.
  auto compute_identifier_view(struct seqinfo_s const & a_sequence) -> View<char> {
    auto const headerlen_signed = static_cast<int>(a_sequence.header_view.size());
    if ((a_sequence.abundance_start == 0) and
        (a_sequence.abundance_end == headerlen_signed)) {
      fatal("Empty sequence identifier.");
    }

    int id_start {0};
    int id_len {0};
    if (a_sequence.abundance_start > 0) {
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


  // Insert id_view into a flat open-addressing dedup table; abort
  // with the offending identifier in the error message if it was
  // already present. Linear probing on collision; the empty default-
  // constructed View<char> marks unused slots (a real identifier
  // cannot be empty - the caller has already aborted on those).
  auto register_unique_identifier(std::vector<View<char>> & hdr_table,
                                  View<char> const id_view) -> void {
    GenericHash<fnv1a> const hdr_hasher;
    auto const hdr_table_size = hdr_table.size();
    auto hdr_idx = hdr_hasher(id_view) % hdr_table_size;
    while (not hdr_table[hdr_idx].empty()) {
      if (hdr_table[hdr_idx] == id_view) {
        std::string const id_str {id_view.data(), id_view.size()};
        fatal("Duplicated sequence identifier: ", id_str);
      }
      hdr_idx = (hdr_idx + 1) % hdr_table_size;
    }
    hdr_table[hdr_idx] = id_view;
  }


  // Open-addressed lookup over seqhashtable, used when d > 0
  // (d = 0 is the dereplication mode and accepts duplicates).
  // Returns true if an identical sequence was already inserted;
  // otherwise records a_sequence at the probed slot and returns false.
  auto is_duplicate_sequence(std::vector<struct seqinfo_s *> & seqhashtable,
                             uint64_t const seqhashsize,
                             struct seqinfo_s & a_sequence) -> bool {
    uint64_t seqhashindex = a_sequence.seqhash % seqhashsize;
    struct seqinfo_s const * seqfound {nullptr};

    while ((seqfound = seqhashtable[seqhashindex]) != nullptr) {
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


  // Pass 1a: header-side population for every entry. Sets up the
  // header_view and (seq, seqlen) pair from the parser entry, then
  // extracts the abundance annotation. Empty identifiers are caught
  // later, in detect_duplicate_identifiers().
  auto index_headers(struct Parameters const & parameters,
                     std::vector<char> const & data_v,
                     std::vector<struct Entry> const & entries,
                     struct Seq_stats & seq_stats,
                     std::vector<struct seqinfo_s> & seqindex_v) -> void {
    Progress progress_hdr("Indexing headers:  ", seq_stats.n_sequences, parameters);
    auto entry_it = entries.cbegin();
    for (auto & a_sequence: seqindex_v) {
        populate_views_from_entry(a_sequence, *entry_it, data_v);

        /* get amplicon abundance */
        find_abundance(a_sequence, seq_stats, entry_it->lineno,
                       parameters.opt_usearch_abundance, parameters.opt_append_abundance);

        progress_hdr.update();
        ++entry_it;
      }
    progress_hdr.done();
  }


  // Pass 1b: detect duplicated identifiers. Walks the populated
  // seqindex_v, computes each identifier subview from the abundance
  // range filled by index_headers(), and inserts it
  // into a local open-addressing dedup table. The table is sized at
  // 2 * n_sequences so the load factor stays <= 0.5: with linear
  // probing this keeps the expected probe length around 1.5 slots,
  // matching the original pre-std::unordered_set implementation. The
  // table is local so its storage is released before the sequence
  // passes allocate their own hashtable. Aborts on a duplicated or
  // empty identifier (compute_identifier_view() catches the latter).
  auto detect_duplicate_identifiers(struct Seq_stats const & seq_stats,
                                    std::vector<struct seqinfo_s> const & seqindex_v,
                                    struct Parameters const & parameters) -> void {
    auto const hdr_table_size = uint64_t{2} * seq_stats.n_sequences;
    std::vector<View<char>> hdr_table(hdr_table_size);

    Progress progress_dup("Checking identifiers:", seq_stats.n_sequences, parameters);
    for (auto const & a_sequence: seqindex_v) {
        register_unique_identifier(hdr_table, compute_identifier_view(a_sequence));
        progress_dup.update();
      }
    progress_dup.done();
  }


  // Pass 2a: compute the Zobrist hash for every sequence. Always run,
  // including in dereplication mode (d = 0), since downstream code
  // relies on a populated seqhash field.
  auto compute_sequence_hashes(Zobrist const & zobrist,
                               struct Seq_stats const & seq_stats,
                               std::vector<struct seqinfo_s> & seqindex_v,
                               struct Parameters const & parameters) -> void {
    Progress progress_hash("Indexing sequences:", seq_stats.n_sequences, parameters);
    for (auto & a_sequence: seqindex_v) {
        a_sequence.seqhash = zobrist.hash(a_sequence.seq, a_sequence.seqlen);
        progress_hash.update();
      }
    progress_hash.done();
  }


  // Pass 2b: detect duplicated sequences via the seqhashtable. Only
  // called when d > 0 (d = 0 is the dereplication mode and accepts
  // duplicates). Stops at the first duplicate; the caller invokes
  // abort_if_duplicated_sequences() afterwards.
  auto detect_duplicate_sequences(struct Seq_stats & seq_stats,
                                  std::vector<struct seqinfo_s> & seqindex_v,
                                  struct Parameters const & parameters) -> void {
    auto const seqhashsize = uint64_t{2} * seq_stats.n_sequences;
    std::vector<struct seqinfo_s *> seqhashtable(seqhashsize);

    Progress progress_dup("Checking duplicates:", seq_stats.n_sequences, parameters);
    for (auto & a_sequence: seqindex_v) {
        if (is_duplicate_sequence(seqhashtable, seqhashsize, a_sequence)) {
            seq_stats.has_duplicates = true;
            break;
          }
        progress_dup.update();
      }
    // Skip done() on the duplicate-detection break: the caller will
    // abort_if_duplicated_sequences() and printing 100% would be misleading.
    if (not seq_stats.has_duplicates) { progress_dup.done(); }
  }


  auto build_index(struct Parameters const & parameters,
                   Zobrist const & zobrist,
                   std::vector<char> const & data_v,
                   std::vector<struct Entry> const & entries,
                   struct Seq_stats & seq_stats,
                   std::vector<struct seqinfo_s> & seqindex_v) -> void {
    seqindex_v.resize(seq_stats.n_sequences);

    index_headers(parameters, data_v, entries, seq_stats, seqindex_v);
    detect_duplicate_identifiers(seq_stats, seqindex_v, parameters);
    compute_sequence_hashes(zobrist, seq_stats, seqindex_v, parameters);
    if (parameters.opt_differences > 0) {
      detect_duplicate_sequences(seq_stats, seqindex_v, parameters);
    }

    abort_if_duplicated_sequences(seq_stats);

    abort_if_missing_abundance(seq_stats);
    sort_index_if_need_be(parameters, seqindex_v);
    print_user_report(parameters, seq_stats);
  }

} // end of anonymous namespace


// ----- class Data -----

Data::Data(struct Parameters const & parameters) {
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


auto Data::info(uint64_t const seqno) const -> struct seqinfo_s const & {
  assert(not seqindex_.empty());  // db_read() / Data ctor must run first
  // bound-check is redundant with -D_GLIBCXX_DEBUG (operator[] is
  // already checked under libstdc++ debug mode), kept here so the
  // precondition is enforced in any assert-enabled build
  assert(seqno < seqindex_.size());
  return seqindex_[seqno];
}


auto Data::sequence_view(uint64_t const seqno) const -> Sequence {
  auto const & rec = info(seqno);
  return {View<char>{rec.seq, nt_bytelength(rec.seqlen)}, rec.seqlen};
}


auto Data::sequence_hash(uint64_t const seqno) const -> uint64_t {
  return info(seqno).seqhash;
}


auto Data::header_view(uint64_t const seqno) const -> View<char> {
  return info(seqno).header_view;
}


auto Data::abundance(uint64_t const seqno) const -> uint64_t {
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
auto Data::fprintseq(std::FILE * stream, unsigned int const seqno) const -> void {
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
                     int64_t const opt_append_abundance) const -> void {
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
                                 bool const opt_usearch_abundance) const -> void {
  auto const & seqinfo = info(seqno);
  auto const * hdrstr = seqinfo.header_view.data();
  auto const hdrlen = static_cast<int>(seqinfo.header_view.size());
  auto const abundance_start = seqinfo.abundance_start;
  auto const abundance_end = seqinfo.abundance_end;

  if (abundance_start < abundance_end) {
      /* print start of header */
      std::fprintf(stream, "%.*s", abundance_start, hdrstr);

      if (opt_usearch_abundance) {
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
                                        bool const opt_usearch_abundance) const -> void {
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

