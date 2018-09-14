/*
    SWARM

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

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

#include "immintrin.h"

#define BLOOM_PATTERN_BITS 10
#define BLOOM_PATTERN_COUNT (1 << BLOOM_PATTERN_BITS)
#define BLOOM_PATTERN_MASK (BLOOM_PATTERN_COUNT - 1)

struct bloom_s
{
  unsigned long mask;
  unsigned long * bitmap;
  unsigned long patterns[BLOOM_PATTERN_COUNT];
};

struct bloom_s * bloom_init(unsigned long size);

void bloom_exit(struct bloom_s * bloom);

inline unsigned long * bloom_adr(struct bloom_s * bloom, unsigned long hash)
{
  return bloom->bitmap + ((hash >> BLOOM_PATTERN_BITS) & bloom->mask);
}

inline void bloom_set(struct bloom_s * bloom, unsigned long hash)
{
  unsigned long p = bloom->patterns[hash & BLOOM_PATTERN_MASK];
  * bloom_adr(bloom, hash) &= ~ p;
}

inline bool bloom_get(struct bloom_s * bloom, unsigned long hash)
{
  unsigned long p = bloom->patterns[hash & BLOOM_PATTERN_MASK];
  return ! (* bloom_adr(bloom, hash) & p);
}
