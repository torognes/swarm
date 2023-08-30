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

#include "utils/nt_codec.h"
#include "variants.h"
#include "zobrist.h"
#include <cstdint>  // uint64_t
#include <cstring>  // std::memcpy


inline void nt_set(char * seq, unsigned int pos, unsigned int base)
{
  // base = replacement nucleotide = encoded as 0, 1, 2, 3
  static constexpr unsigned int divider {5};
  static constexpr unsigned int max_range {31};
  static constexpr unsigned long long int two_bits {3};  // '... 0011' in binary
  const unsigned int whichlong = pos >> divider;
  const uint64_t shift = (pos & max_range) << 1U;  // 0, 2, 4, 6, ..., 60, 62
  const uint64_t mask = two_bits << shift;
  uint64_t mutated_sequence = (reinterpret_cast<uint64_t *>(seq))[whichlong];
  mutated_sequence &= ~ mask;
  mutated_sequence |= (static_cast<uint64_t>(base)) << shift;
  (reinterpret_cast<uint64_t *>(seq))[whichlong] = mutated_sequence;
}

inline void seq_copy(char * seq_a,
                     unsigned int a_start,
                     char * seq_b,
                     unsigned int b_start,
                     unsigned int length)
{
  /* copy part of the compressed sequence b to a */
  for(auto i = 0U; i < length; i++) {
    nt_set(seq_a, a_start + i, nt_extract(seq_b, b_start + i));
  }
}

inline auto seq_identical(char * seq_a,
                          unsigned int a_start,
                          char * seq_b,
                          unsigned int b_start,
                          unsigned int length) -> bool
{
  /* compare parts of two compressed sequences a and b */
  /* return false if different, true if identical */

  bool equal {true};
  for(auto i = 0U; i < length; i++) {
    equal = equal && (nt_extract(seq_a, a_start + i) ==
                      nt_extract(seq_b, b_start + i));
  }
  return equal;
}

void generate_variant_sequence(char * seed_sequence,
                               unsigned int seed_seqlen,
                               struct var_s * var,
                               char * seq,
                               unsigned int * seqlen)
{
  /* generate the actual sequence of a variant */

  switch (var->type)
    {
    case Variant_type::substitution:
      std::memcpy(seq, seed_sequence, nt_bytelength(seed_seqlen));
      nt_set(seq, var->pos, var->base);
      * seqlen = seed_seqlen;
      break;

    case Variant_type::deletion:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos);
      seq_copy(seq, var->pos,
               seed_sequence, var->pos + 1,
               seed_seqlen - var->pos - 1);
      * seqlen = seed_seqlen - 1;
      break;

    case Variant_type::insertion:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos);
      nt_set(seq, var->pos, var->base);
      seq_copy(seq, var->pos + 1,
               seed_sequence, var->pos,
               seed_seqlen - var->pos);
      * seqlen = seed_seqlen + 1;
      break;
    }
}


auto check_variant(char * seed_sequence,
                   unsigned int seed_seqlen,
                   var_s * var,
                   char * amp_sequence,
                   unsigned int amp_seqlen) -> bool
{
  /* make sure seed with given variant is really identical to amp */
  /* we know the hashes are identical */

  bool equal {false};

  switch (var->type)
    {
    case Variant_type::substitution:
      equal = ((seed_seqlen == amp_seqlen) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos)) &&
               (nt_extract(amp_sequence, var->pos) == var->base) &&
               (seq_identical(seed_sequence, var->pos + 1,
                              amp_sequence,  var->pos + 1,
                              seed_seqlen - var->pos - 1)));
      break;

    case Variant_type::deletion:
      equal = (((seed_seqlen - 1) == amp_seqlen) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos)) &&
               (seq_identical(seed_sequence, var->pos + 1,
                              amp_sequence,  var->pos,
                              seed_seqlen - var->pos - 1)));
      break;

    case Variant_type::insertion:
      equal = (((seed_seqlen + 1) == amp_seqlen) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos)) &&
               (nt_extract(amp_sequence, var->pos) == var->base) &&
               (seq_identical(seed_sequence, var->pos,
                              amp_sequence,  var->pos + 1,
                              seed_seqlen - var->pos)));
      break;
    }

  return equal;
}

inline void add_variant(uint64_t hash,
                        Variant_type type,
                        unsigned int pos,
                        unsigned char base,
                        var_s * variant_list,
                        unsigned int * variant_count)
{
#ifdef HASHSTATS
  tries++;
#endif
  var_s * variant = variant_list + (*variant_count)++;
  variant->hash = hash;
  variant->type = type;
  variant->pos = pos;
  variant->base = base;
}

void generate_variants(char * sequence,
                       unsigned int seqlen,
                       uint64_t hash,
                       var_s * variant_list,
                       unsigned int * variant_count)
{
  /* substitutions */

  for(auto offset = 0U; offset < seqlen; offset++)
    {
      const unsigned char current_base = nt_extract(sequence, offset);
      const uint64_t hash1 = hash ^ zobrist_value(offset, current_base);
      for(unsigned char base = 0; base < 4; base++) {
        if (base == current_base) {
          continue;
        }

        const uint64_t hash2 = hash1 ^ zobrist_value(offset, base);
        add_variant(hash2, Variant_type::substitution, offset, base,
                    variant_list, variant_count);

      }
    }

  /* deletions */

  hash = zobrist_hash_delete_first(reinterpret_cast<unsigned char *>(sequence), seqlen);
  add_variant(hash, Variant_type::deletion, 0, 0, variant_list, variant_count);
  unsigned char previous_base = nt_extract(sequence, 0);
  for(auto offset = 1U; offset < seqlen; offset++)
    {
      const unsigned char current_base = nt_extract(sequence, offset);
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
  for(unsigned char base = 0; base < 4; base++)
    {
      const uint64_t hash1 = hash ^ zobrist_value(0, base);
      add_variant(hash1, Variant_type::insertion, 0, base, variant_list, variant_count);
    }
  // insert after each position in the sequence
  for(auto offset = 0U; offset < seqlen; offset++)
    {
      const unsigned char current_base = nt_extract(sequence, offset);
      hash ^= zobrist_value(offset, current_base) ^ zobrist_value(offset + 1, current_base);
      for(unsigned char base = 0; base < 4; base++) {
        if (base == current_base) {
          continue;
        }
        const uint64_t hash1 = hash ^ zobrist_value(offset + 1, base);
        add_variant(hash1, Variant_type::insertion, offset + 1, base,
                    variant_list, variant_count);
      }
    }
}
