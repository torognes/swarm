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

#ifndef SWARM_DB_H
#define SWARM_DB_H

#include "utils/nt_codec.hpp"  // nt_byte_index, nt_extract
#include "utils/seqinfo.hpp"
#include "utils/view.hpp"
#include "utils/zobrist.hpp"
#include <cassert>
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t
#include <memory>  // std::unique_ptr
#include <vector>


struct Parameters;  // defined in swarm.hpp


// Non-owning view of a packed-nucleotide amplicon:
// - length is the nucleotide count
// - encoded views the storage words, with encoded.size() ==
//   nt_wordlength(length) (32 nt per 64-bit word)
// Each sequence starts on a word boundary of the word vector, so
// encoded.data() is 8-byte aligned. Word-wide consumers (seq_identical
// in variants.cpp) index the words directly; byte-level consumers (the
// Zobrist byte-rate zip, Sequence_printer, nucleotide_at below) view
// the words' object representation through encoded.as_bytes() -- the
// same bytes the old byte-typed storage held, since every word arrived
// there via std::memcpy of the same accumulator.
struct Sequence {
  View<uint64_t> encoded;
  unsigned int length;
};


// The nucleotide at 'position' (a nucleotide index, zero-based).
//
// This is the checked way to read a packed sequence, and the check is
// stronger than what a caller can express on its own: 'position' is
// verified against the nucleotide count, where nt_extract() and the
// hand-written bounds it used to need could only see the packed byte
// count -- which also admits the padding nucleotides inside the last
// byte. The View subscript then re-checks the byte index.
inline auto nucleotide_at(Sequence const & sequence, uint64_t const position) -> unsigned char {
  assert(position < sequence.length);
  return nt_extract(sequence.encoded.as_bytes()[nt_byte_index(position)], position);
}


class Data {
public:
  explicit Data(struct Parameters const & parameters);

  auto sequence_count()   const -> unsigned int { return static_cast<unsigned int>(seqindex_.size()); }
  auto longest_sequence() const -> unsigned int { return longest_; }

  auto zobrist() const -> Zobrist const & { return *zobrist_p_; }

  auto info(uint64_t seqno) const -> struct seqinfo_s const &;

  auto sequence_view(uint64_t seqno)   const -> Sequence;
  auto sequence_hash(uint64_t seqno)   const -> uint64_t;
  auto header_view(uint64_t seqno)     const -> View<char>;
  auto abundance(uint64_t seqno)       const -> uint64_t;

private:
  std::vector<char>             data_header_;    // '\0'-terminated headers
  std::vector<uint64_t>         data_sequence_;  // packed sequences, whole 64-bit words
  std::vector<struct seqinfo_s> seqindex_;
  std::unique_ptr<Zobrist>      zobrist_p_;  // deferred: needs longest_sequence
  unsigned int                  longest_ {0};
};


// Writing amplicon labels and sequences used to be four Data member
// functions. They are free now, because none of them ever needed anything
// private: the three label printers read only the four seqinfo_s fields
// that describe the abundance annotation, which info() already hands out,
// and the sequence printer needed a scratch buffer that it can own itself.
// Data is left owning and indexing the database, and no longer writes.
//
// They stay in this translation unit rather than moving to utils/, because
// the annotation format is exactly what the parser here already knows: the
// parser writes abundance_start and abundance_end, and these three read
// them back to strip or replace the annotation. Splitting the two halves
// of one format across two files would cost more than it buys.


// The label as it appeared in the input, with the abundance annotation
// appended when the input carried none and the caller asked for one.
auto fprint_id(std::FILE * stream,
               struct seqinfo_s const & seqinfo,
               bool opt_usearch_abundance,
               int64_t opt_append_abundance) -> void;

// The label with its abundance annotation removed.
auto fprint_id_noabundance(std::FILE * stream,
                           struct seqinfo_s const & seqinfo,
                           bool opt_usearch_abundance) -> void;

// The label with its abundance annotation replaced by 'new_abundance'.
auto fprint_id_with_new_abundance(std::FILE * stream,
                                  struct seqinfo_s const & seqinfo,
                                  uint64_t new_abundance,
                                  bool opt_usearch_abundance) -> void;


// Decodes packed nucleotides to ascii and writes them, one FASTA sequence
// line per call.
//
// A class rather than a free function because the decode needs scratch
// space, and that space needs a size rule the callers should not have to
// know: a packed byte decodes to a whole group of four characters, so the
// last group of a sequence whose length is not a multiple of four runs up
// to three characters past it. Sizing the buffer from longest_sequence()
// *without* rounding up wrote out of bounds on the very first -w run.
// Owning the buffer keeps that rule next to the code that depends on it,
// and replaces a mutable member of Data -- a const member function that
// quietly wrote to its own object.
//
// Construct once per output file and reuse, as the seeds writers do with
// NwAligner; the buffer is sized for the longest sequence in the database.
class Sequence_printer {
public:
  explicit Sequence_printer(unsigned int longest_sequence);

  // Named print(), not fprint(): a member named fprint would hide the free
  // fprint() overloads inside this class, and the body below calls one.
  auto print(std::FILE * stream, Sequence const & sequence) const -> void;

private:
  mutable std::vector<char> decode_buffer_;
};

#endif  // SWARM_DB_H


