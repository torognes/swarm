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

/*
  Blocked bloom filter with precomputed bit patterns
  as described in

  Putze F, Sanders P, Singler J (2009)
  Cache-, Hash- and Space-Efficient Bloom Filters
  Journal of Experimental Algorithmics, 14, 4
  https://doi.org/10.1145/1498698.1594230
*/

#include "swarm.h"

void bloom_patterns_generate(struct bloom_s * b)
{
  const unsigned int k = 8;
  for (unsigned int i = 0; i < BLOOM_PATTERN_COUNT; i++)
    {
      unsigned long pattern = 0;
      for (unsigned int j = 0; j < k; j++)
        {
          unsigned long onebit;
          onebit = 1ULL << (random() & 63);
          while (pattern & onebit)
            onebit = 1ULL << (random() & 63);
          pattern |= onebit;
        }
      b->patterns[i] = pattern;
    }
}

void bloom_zap(struct bloom_s * b)
{
  memset(b->bitmap, 0xff, b->size);
}

struct bloom_s * bloom_init(unsigned long size)
{
  // Size is in bytes for full bitmap, must be power of 2
  // at least 8
  size = MAX(size, 8);

  struct bloom_s * b = (struct bloom_s *) xmalloc(sizeof(struct bloom_s));

  b->size = size;

  b->mask = (size >> 3) - 1;

  b->bitmap = (unsigned long *) xmalloc(size);

  bloom_zap(b);

  bloom_patterns_generate(b);

  return b;
}

void bloom_exit(struct bloom_s * b)
{
  free(b->bitmap);
  free(b);
}
