/*
    SWARM

    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

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

/* Variant information */

constexpr int substitution {1};
constexpr int deletion {2};
constexpr int insertion {3};

struct var_s
{
  uint64_t hash;
  unsigned int pos;
  unsigned char type;
  unsigned char base;
  unsigned short dummy; /* for alignment padding only */
};

void generate_variant_sequence(char * seed_sequence,
                               unsigned int seed_seqlen,
                               struct var_s * var,
                               char * seq,
                               unsigned int * seqlen);

auto check_variant(char * seed_sequence,
                   unsigned int seed_seqlen,
                   struct var_s * var,
                   char * amp_sequence,
                   unsigned int amp_seqlen) -> bool;

void generate_variants(char * sequence,
                       unsigned int seqlen,
                       uint64_t hash,
                       struct var_s * variant_list,
                       unsigned int * variant_count);
