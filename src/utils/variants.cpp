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

#include "../db.hpp"
#include "nt_codec.hpp"
#include "span.hpp"
#include "view.hpp"
#include "variants.hpp"
#include <algorithm>  // std::copy
#include <cassert>  // assert
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // uint64_t
#include <cstring>  // std::memcpy
#include <iterator>  // std::next
#include <vector>



namespace {

#ifndef NDEBUG
  // C++17 refactoring: [[maybe_unused]]
  constexpr std::size_t nt_per_byte = 4;  // 4 nucleotides packed per byte
#endif


  inline auto nt_set(Span<char> const seq, unsigned int const pos, unsigned int const base) -> void
  {
    // base = replacement nucleotide = encoded as 0, 1, 2, 3
    static constexpr auto divider = 5U;
    static constexpr auto max_range = 31U;
    static constexpr auto two_bits = 3ULL;  // '... 0011' in binary
    const auto whichlong = pos >> divider;
    const uint64_t shift = static_cast<uint64_t>(pos & max_range) << 1U;  // 0, 2, 4, 6, ..., 60, 62
    const uint64_t mask = compl (two_bits << shift);
    // read-modify-write the target 64-bit word. std::memcpy avoids the
    // strict-aliasing undefined behaviour of punning a char buffer through
    // a uint64_t* (see C++ Weekly #185); optimisers fold the round-trip
    // back into a single load and store.
    // C++20 refactoring: std::bit_cast
    const auto byte_offset = static_cast<std::ptrdiff_t>(whichlong)
                           * static_cast<std::ptrdiff_t>(sizeof(uint64_t));
    char * const target_word = std::next(seq.data(), byte_offset);
    uint64_t mutated_position {0};
    std::memcpy(&mutated_position, target_word, sizeof(mutated_position));
    mutated_position &= mask;
    mutated_position |= (static_cast<uint64_t>(base)) << shift;
    std::memcpy(target_word, &mutated_position, sizeof(mutated_position));
  }


  inline auto seq_copy(Span<char> const seq_a,
                       unsigned int a_start,
                       View<char> seq_b,
                       unsigned int b_start,
                       unsigned int length) -> void
  {
    /* copy part of the compressed sequence b to a */
    //
    // One read-modify-write of a whole word per nucleotide, which looks
    // like the write-side counterpart of seq_identical below and is not:
    // writing a destination word per iteration instead, with nt_window()
    // for the source and a mask for the partial first and last words,
    // was built and measured 4.5 % slower at -d 1 -f (18SV9 10%, user
    // CPU, single-threaded, five alternating pairs: +3.29, +1.91, +3.02,
    // +3.02, +1.15 s on a 53.4 s baseline, both orders agreeing). It
    // does thirty times less work in the loop and still loses: inlined
    // into the fastidious worker it added 664 bytes to it, the same way
    // the packed-byte block cost search8 and search16 (see fill_channel
    // in dseq_fill.hpp). seq_identical is where that machinery pays.
    assert(static_cast<std::size_t>(a_start) + length <= seq_a.size() * nt_per_byte);
    assert(static_cast<std::size_t>(b_start) + length <= seq_b.size() * nt_per_byte);
    for (auto i = 0U; i < length; ++i) {
      nt_set(seq_a, a_start + i, nt_extract(seq_b[nt_byte_index(b_start + i)], b_start + i));
    }
  }


  constexpr unsigned int bits_per_nt = 2;   // 2-bit fields
  constexpr unsigned int nt_per_word = 32;  // 32 nucleotides in a uint64_t
  constexpr unsigned int bits_per_word = bits_per_nt * nt_per_word;


  // The 64-bit word at 'word_index' of a packed sequence. Every packed
  // sequence occupies a whole number of these: nt_bytelength() rounds
  // the byte count up to a multiple of 8, so a word read stays inside
  // the view for any nucleotide position the view covers.
  //
  // std::memcpy avoids the strict-aliasing undefined behaviour of
  // punning a char buffer through a uint64_t* (see nt_set above);
  // optimisers fold it back into a single load.
  // C++20 refactoring: std::bit_cast
  inline auto packed_word(View<char> const seq, uint64_t const word_index) -> uint64_t
  {
    assert(seq.size() % sizeof(uint64_t) == 0);
    assert((word_index + 1) * sizeof(uint64_t) <= seq.size());
    auto const byte_offset = static_cast<std::ptrdiff_t>(word_index * sizeof(uint64_t));
    uint64_t word {0};
    std::memcpy(&word, std::next(seq.data(), byte_offset), sizeof(word));
    return word;
  }


  // The 32 nucleotides starting at 'position', packed lowest field
  // first. A variant sits anywhere in the sequence, so 'position' is
  // not word-aligned in general and the window is assembled from the
  // word holding it and, when it straddles two, the one after.
  inline auto nt_window(View<char> const seq, uint64_t const position) -> uint64_t
  {
    auto const word_index = position / nt_per_word;
    auto const shift = bits_per_nt * (position % nt_per_word);
    auto const word = packed_word(seq, word_index) >> shift;

    // an aligned window is already complete, and shifting a 64-bit
    // value by 64 is undefined; a window opening in the last word has
    // no successor to draw its high fields from, and does not need one
    // (they lie past the end of the sequence, so no caller compares
    // them -- see the mask in seq_identical)
    auto const has_next_word = ((word_index + 2) * sizeof(uint64_t)) <= seq.size();
    if ((shift == 0) or (not has_next_word)) {
      return word;
    }
    return word | (packed_word(seq, word_index + 1) << (bits_per_word - shift));
  }


  inline auto seq_identical(View<char> seq_a,
                            unsigned int a_start,
                            View<char> seq_b,
                            unsigned int b_start,
                            unsigned int length) -> bool
  {
    /* compare parts of two compressed sequences a and b */
    /* return false if different, true if identical */
    assert(static_cast<std::size_t>(a_start) + length <= seq_a.size() * nt_per_byte);
    assert(static_cast<std::size_t>(b_start) + length <= seq_b.size() * nt_per_byte);

    // 32 nucleotides per iteration. The two windows carry the same
    // nucleotides at the same field positions whatever a_start and
    // b_start are -- equal for a substitution, one apart for an
    // insertion or a deletion -- so one xor compares all 32.
    for (auto compared = 0U; compared < length; compared += nt_per_word) {
      auto const difference = nt_window(seq_a, a_start + compared)
                            ^ nt_window(seq_b, b_start + compared);

      // The final window is partial. Its fields past 'length' are
      // padding: zero in a parsed sequence, but stale in the variant
      // buffer, which generate_variant_sequence rewrites in place and
      // which the fastidious path passes here. They are masked out
      // rather than compared.
      auto const remaining = length - compared;
      if (remaining < nt_per_word) {
        auto const mask = (uint64_t {1} << (bits_per_nt * remaining)) - 1;
        return (difference & mask) == 0;
      }
      if (difference != 0) {
        return false;
      }
    }
    return true;
  }


  inline auto add_variant(uint64_t hash,
                          Variant_type type,
                          unsigned int pos,
                          unsigned char base,
                          std::vector<struct var_s>& variant_list,
                          unsigned int & variant_count) -> void
  {
    var_s & variant = variant_list[variant_count];
    ++variant_count;
    variant.hash = hash;
    variant.type = type;
    variant.pos = pos;
    variant.base = base;
  }

}  // anonymous namespace


auto generate_variant_sequence(Sequence const & seed,
                               struct var_s const & var,
                               std::vector<char> & buffer) -> Sequence
{
  /* generate the actual sequence of a variant */

  auto const seed_seqlen = seed.length;
  auto const seq_span = make_span(buffer);
  auto seqlen = 0U;

  switch (var.type)
    {
    case Variant_type::substitution:
      std::copy(seed.encoded.cbegin(), seed.encoded.cend(), buffer.begin());
      nt_set(seq_span, var.pos, var.base);
      seqlen = seed_seqlen;
      break;

    case Variant_type::deletion:
      seq_copy(seq_span, 0,
               seed.encoded, 0,
               var.pos);
      seq_copy(seq_span, var.pos,
               seed.encoded, var.pos + 1,
               seed_seqlen - var.pos - 1);
      seqlen = seed_seqlen - 1;
      break;

    case Variant_type::insertion:
      seq_copy(seq_span, 0,
               seed.encoded, 0,
               var.pos);
      nt_set(seq_span, var.pos, var.base);
      seq_copy(seq_span, var.pos + 1,
               seed.encoded, var.pos,
               seed_seqlen - var.pos);
      seqlen = seed_seqlen + 1;
      break;
    }

  // a view of the caller's buffer, not of storage this function owns:
  // see the note on the declaration
  assert(seqlen != 0);
  return Sequence{make_view(buffer).first(nt_bytelength(seqlen)), seqlen};
}


auto check_variant(Sequence const & seed,
                   struct var_s const & var,
                   Sequence const & amp) -> bool
{
  /* make sure seed with given variant is really identical to amp */
  /* we know the hashes are identical */

  auto const seed_seqlen = seed.length;
  auto const amp_seqlen = amp.length;

  bool equal {false};

  switch (var.type)
    {
    case Variant_type::substitution:
      equal = ((seed_seqlen == amp_seqlen) and
               seq_identical(seed.encoded, 0,
                             amp.encoded, 0,
                             var.pos) and
               (nucleotide_at(amp, var.pos) == var.base) and
               seq_identical(seed.encoded, var.pos + 1,
                             amp.encoded,  var.pos + 1,
                             seed_seqlen - var.pos - 1));
      break;

    case Variant_type::deletion:
      equal = (((seed_seqlen - 1) == amp_seqlen) and
               seq_identical(seed.encoded, 0,
                             amp.encoded, 0,
                             var.pos) and
               seq_identical(seed.encoded, var.pos + 1,
                             amp.encoded,  var.pos,
                             seed_seqlen - var.pos - 1));
      break;

    case Variant_type::insertion:
      equal = (((seed_seqlen + 1) == amp_seqlen) and
               seq_identical(seed.encoded, 0,
                             amp.encoded, 0,
                             var.pos) and
               (nucleotide_at(amp, var.pos) == var.base) and
               seq_identical(seed.encoded, var.pos,
                             amp.encoded,  var.pos + 1,
                             seed_seqlen - var.pos));
      break;
    }

  return equal;
}


auto generate_variants(Zobrist const & zobrist,
                       Sequence const & seq,
                       uint64_t hash,
                       std::vector<struct var_s>& variant_list) -> View<struct var_s>
{
  auto const seqlen = seq.length;

  auto variant_count = 0U;
  /* substitutions */

  for (auto position = 0U; position < seqlen; ++position)
    {
      const auto current_base = nucleotide_at(seq, position);
      const auto hash1 = hash ^ zobrist.value(position, current_base);
      for (unsigned char base = 0; base < 4; ++base) {
        if (base == current_base) {
          continue;
        }

        const auto hash2 = hash1 ^ zobrist.value(position, base);
        add_variant(hash2, Variant_type::substitution, position, base,
                    variant_list, variant_count);

      }
    }

  /* deletions */

  hash = zobrist.hash_delete_first(seq);
  add_variant(hash, Variant_type::deletion, 0, 0, variant_list, variant_count);
  auto previous_base = nucleotide_at(seq, 0);
  for (auto offset = 1U; offset < seqlen; ++offset)
    {
      const auto current_base = nucleotide_at(seq, offset);
      if (current_base == previous_base) {
        continue;
      }
      hash ^= zobrist.value(offset - 1, previous_base) ^ zobrist.value(offset - 1, current_base);
      add_variant(hash, Variant_type::deletion, offset, 0, variant_list, variant_count);
      previous_base = current_base;
    }

  /* insertions */

  hash = zobrist.hash_insert_first(seq);
  // insert before the first position in the sequence
  for (unsigned char base = 0; base < 4; ++base)
    {
      const auto hash1 = hash ^ zobrist.value(0, base);
      add_variant(hash1, Variant_type::insertion, 0, base, variant_list, variant_count);
    }
  // insert after each position in the sequence
  for (auto position = 0U; position < seqlen; ++position)
    {
      const auto current_base = nucleotide_at(seq, position);
      hash ^= zobrist.value(position, current_base) ^ zobrist.value(position + 1, current_base);
      for (unsigned char base = 0; base < 4; ++base) {
        if (base == current_base) {
          continue;
        }
        const auto hash1 = hash ^ zobrist.value(position + 1, base);
        add_variant(hash1, Variant_type::insertion, position + 1, base,
                    variant_list, variant_count);
      }
    }

  return make_view(variant_list).first(variant_count);
}
