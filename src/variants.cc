/*
    SWARM

    Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe

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

inline void nt_set(char * seq, unsigned int pos, unsigned int base)
{
  unsigned int whichlong = pos >> 5;
  uint64_t shift = (pos & 31) << 1;
  uint64_t mask = 3ULL << shift;
  uint64_t x = ((uint64_t *) seq)[whichlong];
  x &= ~ mask;
  x |= ((uint64_t)base) << shift;
  ((uint64_t *) seq)[whichlong] = x;
}

inline void seq_copy(char * a,
                     int a_start,
                     char * b,
                     int b_start,
                     int length)
{
  /* copy part of the compressed sequence b to a */
  for(int i = 0; i < length; i++)
    nt_set(a, a_start + i, nt_extract(b, b_start + i));
}

inline bool seq_identical(char * a,
                          int a_start,
                          char * b,
                          int b_start,
                          int length)
{
  /* compare parts of two compressed sequences a and b */
  /* return false if different, true if identical */

  for(int i = 0; i < length; i++)
    if (nt_extract(a, a_start + i) != nt_extract(b, b_start + i))
      return false;
  return true;
}

void generate_variant_sequence(char * seed_sequence,
                               unsigned int seed_seqlen,
                               struct var_s * var,
                               char * seq,
                               int * seqlen)
{
  /* generate the actual sequence of a variant */

  switch (var->type)
    {
    case identical:
      memcpy(seq, seed_sequence, nt_bytelength(seed_seqlen));
      * seqlen = seed_seqlen;
      break;

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

    default:
      fatal("Unknown variant");
    }
}


bool check_variant(char * seed_sequence,
                   unsigned int seed_seqlen,
                   var_s * var,
                   char * amp_sequence,
                   unsigned int amp_seqlen)
{
  /* make sure seed with given variant is really identical to amp */
  /* we know the hashes are identical */

  switch (var->type)
    {
    case identical:
      if (seed_seqlen != amp_seqlen)
        return false;
      return seq_identical(seed_sequence, 0,
                           amp_sequence, 0,
                           seed_seqlen);

    case substitution:
      if (seed_seqlen != amp_seqlen)
        return false;
      if (! seq_identical(seed_sequence, 0,
                          amp_sequence, 0,
                          var->pos))
        return false;
      if (nt_extract(amp_sequence, var->pos) != var->base)
        return false;
      return seq_identical(seed_sequence, var->pos + 1,
                           amp_sequence,  var->pos + 1,
                           seed_seqlen - var->pos - 1);

    case deletion:
      if ((seed_seqlen - 1) != amp_seqlen)
        return false;
      if (! seq_identical(seed_sequence, 0,
                          amp_sequence, 0,
                          var->pos))
        return false;
      return seq_identical(seed_sequence, var->pos + 1,
                           amp_sequence,  var->pos,
                           seed_seqlen - var->pos - 1);

    case insertion:
      if ((seed_seqlen + 1) != amp_seqlen)
        return false;
      if (! seq_identical(seed_sequence, 0,
                          amp_sequence, 0,
                          var->pos))
        return false;
      if (nt_extract(amp_sequence, var->pos) != var->base)
        return false;
      return seq_identical(seed_sequence, var->pos,
                           amp_sequence,  var->pos + 1,
                           seed_seqlen - var->pos);

    default:
      fatal("Unknown variant");
      return false;
    }
}

inline void add_variant(uint64_t hash,
                        unsigned char type,
                        unsigned int pos,
                        unsigned int base,
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
                       unsigned int * variant_count,
                       bool include_identical)
{
  /* identical non-variant */

  if (include_identical)
    add_variant(hash, identical, 0, 0, variant_list, variant_count);

  /* substitutions */

  for(unsigned int i = 0; i < seqlen; i++)
    {
      unsigned int base = nt_extract(sequence, i);
      uint64_t hash1 = hash ^ zobrist_value(i, base);
      for (unsigned int v = 0; v < 4; v ++)
        if (v != base)
          {
            uint64_t hash2 = hash1 ^ zobrist_value(i, v);
            add_variant(hash2, substitution, i, v,
                        variant_list, variant_count);
          }
    }

  /* deletions */

  hash = zobrist_hash_delete_first((unsigned char *) sequence, seqlen);
  add_variant(hash, deletion, 0, 0, variant_list, variant_count);
  unsigned int base_deleted = nt_extract(sequence, 0);
  for(unsigned int i = 1; i < seqlen; i++)
    {
      unsigned int v = nt_extract(sequence, i);
      if (v != base_deleted)
        {
          hash ^= zobrist_value(i - 1, base_deleted) ^ zobrist_value(i - 1, v);
          add_variant(hash, deletion, i, 0, variant_list, variant_count);
          base_deleted = v;
        }
    }

  /* insertions */

  hash = zobrist_hash_insert_first((unsigned char *) sequence, seqlen);
  for (unsigned int v = 0; v < 4; v++)
    {
      uint64_t hash1 = hash ^ zobrist_value(0, v);
      add_variant(hash1, insertion, 0, v, variant_list, variant_count);
    }
  for (unsigned int i = 0; i < seqlen; i++)
    {
      unsigned int base = nt_extract(sequence, i);
      hash ^= zobrist_value(i, base) ^ zobrist_value(i+1, base);
      for (unsigned int v = 0; v < 4; v++)
        if (v != base)
          {
            uint64_t hash1 = hash ^ zobrist_value(i + 1, v);
            add_variant(hash1, insertion, i + 1, v,
                        variant_list, variant_count);
          }
    }
}
