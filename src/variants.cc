/*
    SWARM

    Copyright (C) 2012-2025 Torbjorn Rognes and Frederic Mahe

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

#include "utils/nt_codec.h"
#include "variants.h"
#include "zobrist.h"
#include <cstdint>  // uint64_t
#include <cstring>  // std::memcpy
#include <iterator>  // std::next
#include <vector>


inline auto nt_set(char * const seq, unsigned int const pos, unsigned int const base) -> void
{
  // base = replacement nucleotide = encoded as 0, 1, 2, 3
  static constexpr auto divider = 5U;
  static constexpr auto max_range = 31U;
  static constexpr auto two_bits = 3ULL;  // '... 0011' in binary
  const auto whichlong = pos >> divider;
  const uint64_t shift = (pos & max_range) << 1U;  // 0, 2, 4, 6, ..., 60, 62
  const uint64_t mask = compl (two_bits << shift);
  auto & mutated_position = *std::next(reinterpret_cast<uint64_t *>(seq), whichlong);
  mutated_position &= mask;
  mutated_position |= (static_cast<uint64_t>(base)) << shift;
}


inline auto seq_copy(char * seq_a,
                     unsigned int a_start,
                     char const * seq_b,
                     unsigned int b_start,
                     unsigned int length) -> void
{
  /* copy part of the compressed sequence b to a */
  for(auto i = 0U; i < length; ++i) {
    nt_set(seq_a, a_start + i, nt_extract(seq_b, b_start + i));
  }
}


inline auto seq_identical(char const * seq_a,
                          unsigned int a_start,
                          char const * seq_b,
                          unsigned int b_start,
                          unsigned int length) -> bool
{
  /* compare parts of two compressed sequences a and b */
  /* return false if different, true if identical */
  for(auto i = 0U; i < length; ++i) {
    if (nt_extract(seq_a, a_start + i) != nt_extract(seq_b, b_start + i)) {
      return false;
    }
  }
  return true;
}


auto generate_variant_sequence(char * seed_sequence,
                               unsigned int seed_seqlen,
                               struct var_s & var,
                               std::vector<char>& seq,
                               unsigned int & seqlen) -> void
{
  /* generate the actual sequence of a variant */

  switch (var.type)
    {
    case Variant_type::substitution:
      std::memcpy(seq.data(), seed_sequence, nt_bytelength(seed_seqlen));
      nt_set(seq.data(), var.pos, var.base);
      seqlen = seed_seqlen;
      break;

    case Variant_type::deletion:
      seq_copy(seq.data(), 0,
               seed_sequence, 0,
               var.pos);
      seq_copy(seq.data(), var.pos,
               seed_sequence, var.pos + 1,
               seed_seqlen - var.pos - 1);
      seqlen = seed_seqlen - 1;
      break;

    case Variant_type::insertion:
      seq_copy(seq.data(), 0,
               seed_sequence, 0,
               var.pos);
      nt_set(seq.data(), var.pos, var.base);
      seq_copy(seq.data(), var.pos + 1,
               seed_sequence, var.pos,
               seed_seqlen - var.pos);
      seqlen = seed_seqlen + 1;
      break;
    }
}


auto check_variant(char const * seed_sequence,
                   unsigned int seed_seqlen,
                   struct var_s & var,
                   char const * amp_sequence,
                   unsigned int amp_seqlen) -> bool
{
  /* make sure seed with given variant is really identical to amp */
  /* we know the hashes are identical */

  bool equal {false};

  switch (var.type)
    {
    case Variant_type::substitution:
      equal = ((seed_seqlen == amp_seqlen) and
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var.pos)) and
               (nt_extract(amp_sequence, var.pos) == var.base) and
               (seq_identical(seed_sequence, var.pos + 1,
                              amp_sequence,  var.pos + 1,
                              seed_seqlen - var.pos - 1)));
      break;

    case Variant_type::deletion:
      equal = (((seed_seqlen - 1) == amp_seqlen) and
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var.pos)) and
               (seq_identical(seed_sequence, var.pos + 1,
                              amp_sequence,  var.pos,
                              seed_seqlen - var.pos - 1)));
      break;

    case Variant_type::insertion:
      equal = (((seed_seqlen + 1) == amp_seqlen) and
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var.pos)) and
               (nt_extract(amp_sequence, var.pos) == var.base) and
               (seq_identical(seed_sequence, var.pos,
                              amp_sequence,  var.pos + 1,
                              seed_seqlen - var.pos)));
      break;
    }

  return equal;
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


auto generate_variants(char * sequence,
                       unsigned int seqlen,
                       uint64_t hash,
                       std::vector<struct var_s>& variant_list) -> unsigned int
{
  auto variant_count = 0U;
  /* substitutions */

  for(auto position = 0U; position < seqlen; ++position)
    {
      const auto current_base = nt_extract(sequence, position);
      const auto hash1 = hash ^ zobrist_value(position, current_base);
      for(unsigned char base = 0; base < 4; ++base) {
        if (base == current_base) {
          continue;
        }

        const auto hash2 = hash1 ^ zobrist_value(position, base);
        add_variant(hash2, Variant_type::substitution, position, base,
                    variant_list, variant_count);

      }
    }

  /* deletions */

  hash = zobrist_hash_delete_first(sequence, seqlen);
  add_variant(hash, Variant_type::deletion, 0, 0, variant_list, variant_count);
  auto previous_base = nt_extract(sequence, 0);
  for(auto offset = 1U; offset < seqlen; ++offset)
    {
      const auto current_base = nt_extract(sequence, offset);
      if (current_base == previous_base) {
        continue;
      }
      hash ^= zobrist_value(offset - 1, previous_base) ^ zobrist_value(offset - 1, current_base);
      add_variant(hash, Variant_type::deletion, offset, 0, variant_list, variant_count);
      previous_base = current_base;
    }

  /* insertions */

  hash = zobrist_hash_insert_first(reinterpret_cast<unsigned char *>(sequence), seqlen);
  // insert before the first position in the sequence
  for(unsigned char base = 0; base < 4; ++base)
    {
      const auto hash1 = hash ^ zobrist_value(0, base);
      add_variant(hash1, Variant_type::insertion, 0, base, variant_list, variant_count);
    }
  // insert after each position in the sequence
  for(auto position = 0U; position < seqlen; ++position)
    {
      const auto current_base = nt_extract(sequence, position);
      hash ^= zobrist_value(position, current_base) ^ zobrist_value(position + 1, current_base);
      for(unsigned char base = 0; base < 4; ++base) {
        if (base == current_base) {
          continue;
        }
        const auto hash1 = hash ^ zobrist_value(position + 1, base);
        add_variant(hash1, Variant_type::insertion, position + 1, base,
                    variant_list, variant_count);
      }
    }

  return variant_count;
}
