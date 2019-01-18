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

#include "immintrin.h"

struct bloomflex_s
{
  unsigned long size; /* size in number of longs (8 bytes) */
  unsigned long pattern_shift;
  unsigned long pattern_count;
  unsigned long pattern_mask;
  unsigned long pattern_k;
  unsigned long * bitmap;
  unsigned long * patterns;
};

struct bloomflex_s * bloomflex_init(unsigned long size, unsigned int k);

void bloomflex_exit(struct bloomflex_s * b);

inline unsigned long * bloomflex_adr(struct bloomflex_s * b, unsigned long h)
{
  return b->bitmap + ((h >> b->pattern_shift) % b->size);
}

inline unsigned long bloomflex_pat(struct bloomflex_s * b,
                                     unsigned long h)
{
  return b->patterns[h & b->pattern_mask];
}

inline void bloomflex_set(struct bloomflex_s * b, unsigned long h)
{
  * bloomflex_adr(b, h) &= ~ bloomflex_pat(b, h);
}

inline bool bloomflex_get(struct bloomflex_s * b, unsigned long h)
{
  return ! (* bloomflex_adr(b, h) & bloomflex_pat(b, h));
}
