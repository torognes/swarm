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
#include "variants.h"
#include "zobrist.h"

inline void nt_set(char * seq, unsigned int pos, unsigned int base)
{
  static constexpr unsigned int divider {5};
  static constexpr unsigned int max_range {31};
  static constexpr unsigned long long int two_bits {3};  // '... 0011' in binary
  unsigned int whichlong = pos >> divider;
  uint64_t shift = (pos & max_range) << 1;  // 0, 2, 4, 6, ..., 60, 62
  uint64_t mask = two_bits << shift;
  uint64_t x = (reinterpret_cast<uint64_t *>(seq))[whichlong];
  x &= ~ mask;
  x |= (static_cast<uint64_t>(base)) << shift;
  (reinterpret_cast<uint64_t *>(seq))[whichlong] = x;
}

inline void seq_copy(char * a,
                     unsigned int a_start,
                     char * b,
                     unsigned int b_start,
                     unsigned int length)
{
  /* copy part of the compressed sequence b to a */
  for(auto i = 0U; i < length; i++) {
    nt_set(a, a_start + i, nt_extract(b, b_start + i));
  }
}

inline auto seq_identical(char * a,
                          unsigned int a_start,
                          char * b,
                          unsigned int b_start,
                          unsigned int length) -> bool
{
  /* compare parts of two compressed sequences a and b */
  /* return false if different, true if identical */

  bool equal {true};
  for(auto i = 0U; i < length; i++) {
    equal = equal && (nt_extract(a, a_start + i) ==
                      nt_extract(b, b_start + i));
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

  assert((var->type == substitution) ||
         (var->type == deletion) ||
         (var->type == insertion));

  switch (var->type)
    {
    case substitution:
      memcpy(seq, seed_sequence, nt_bytelength(seed_seqlen));
      nt_set(seq, var->pos, var->base);
      * seqlen = seed_seqlen;
      break;

    case deletion:
      seq_copy(seq, 0,
               seed_sequence, 0,
               var->pos);
      seq_copy(seq, var->pos,
               seed_sequence, var->pos + 1,
               seed_seqlen - var->pos - 1);
      * seqlen = seed_seqlen - 1;
      break;

    case insertion:
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

  assert((var->type == substitution) ||
         (var->type == deletion) ||
         (var->type == insertion));

  bool equal {false};

  switch (var->type)
    {
    case substitution:
      equal = ((seed_seqlen == amp_seqlen) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos)) &&
               (nt_extract(amp_sequence, var->pos) == var->base) &&
               (seq_identical(seed_sequence, var->pos + 1,
                              amp_sequence,  var->pos + 1,
                              seed_seqlen - var->pos - 1)));
      break;

    case deletion:
      equal = (((seed_seqlen - 1) == amp_seqlen) &&
               (seq_identical(seed_sequence, 0,
                              amp_sequence, 0,
                              var->pos)) &&
               (seq_identical(seed_sequence, var->pos + 1,
                              amp_sequence,  var->pos,
                              seed_seqlen - var->pos - 1)));
      break;

    case insertion:
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
                        unsigned char type,
                        unsigned int pos,
                        unsigned char base,
                        var_s * variant_list,
                        unsigned int * variant_count)
{
#ifdef HASHSTATS
  tries++;
#endif
  var_s * v = variant_list + (*variant_count)++;
  v->hash = hash;
  v->type = type;
  v->pos = pos;
  v->base = base;
}

void generate_variants(char * sequence,
                       unsigned int seqlen,
                       uint64_t hash,
                       var_s * variant_list,
                       unsigned int * variant_count)
{
  /* substitutions */

  for(auto i = 0U; i < seqlen; i++)
    {
      unsigned char base = nt_extract(sequence, i);
      uint64_t hash1 = hash ^ zobrist_value(i, base);
      for(unsigned char v = 0; v < 4; v ++) {
        if (v != base)
          {
            uint64_t hash2 = hash1 ^ zobrist_value(i, v);
            add_variant(hash2, substitution, i, v,
                        variant_list, variant_count);
          }
      }
    }

  /* deletions */

  hash = zobrist_hash_delete_first(reinterpret_cast<unsigned char *>(sequence), seqlen);
  add_variant(hash, deletion, 0, 0, variant_list, variant_count);
  unsigned char base_deleted = nt_extract(sequence, 0);
  for(auto i = 1U; i < seqlen; i++)
    {
      unsigned char v = nt_extract(sequence, i);
      if (v != base_deleted)
        {
          hash ^= zobrist_value(i - 1, base_deleted) ^ zobrist_value(i - 1, v);
          add_variant(hash, deletion, i, 0, variant_list, variant_count);
          base_deleted = v;
        }
    }

  /* insertions */

  hash = zobrist_hash_insert_first(reinterpret_cast<unsigned char *>(sequence), seqlen);
  for(unsigned char v = 0; v < 4; v++)
    {
      uint64_t hash1 = hash ^ zobrist_value(0, v);
      add_variant(hash1, insertion, 0, v, variant_list, variant_count);
    }
  for(auto i = 0U; i < seqlen; i++)
    {
      unsigned char base = nt_extract(sequence, i);
      hash ^= zobrist_value(i, base) ^ zobrist_value(i+1, base);
      for(unsigned char v = 0; v < 4; v++) {
        if (v != base)
          {
            uint64_t hash1 = hash ^ zobrist_value(i + 1, v);
            add_variant(hash1, insertion, i + 1, v,
                        variant_list, variant_count);
          }
      }
    }
}
