/*
    SWARM

    Copyright (C) 2012-2020 Torbjorn Rognes and Frederic Mahe

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

void bloomflex_patterns_generate(struct bloomflex_s * b);

void bloomflex_patterns_generate(struct bloomflex_s * b)
{
#if 0
  printf("Generating %" PRIu64 " patterns with %" PRIu64 " bits set.\n",
         b->pattern_count,
         b->pattern_k);
#endif
  for(auto i = 0U; i < b->pattern_count; i++)
    {
      uint64_t pattern {0};
      for(auto j = 0U; j < b->pattern_k; j++)
        {
          uint64_t onebit {0};
          onebit = 1ULL << (arch_random() & 63);
          while ((pattern & onebit) != 0U) {
            onebit = 1ULL << (arch_random() & 63);
          }
          pattern |= onebit;
        }
      b->patterns[i] = pattern;
    }
}

struct bloomflex_s * bloomflex_init(uint64_t size, unsigned int k)
{
  /* Input size is in bytes for full bitmap */

  auto * b = static_cast<struct bloomflex_s *>(xmalloc(sizeof(struct bloomflex_s)));
  b->size = size >> 3;

  b->pattern_shift = 16;
  b->pattern_count = 1 << b->pattern_shift;
  b->pattern_mask = b->pattern_count - 1;
  b->pattern_k = k;

  b->patterns = static_cast<uint64_t *>(xmalloc(b->pattern_count * 8));
  bloomflex_patterns_generate(b);

  b->bitmap = static_cast<uint64_t *>(xmalloc(size));
  memset(b->bitmap, 0xff, size);

  return b;
}

void bloomflex_exit(struct bloomflex_s * b)
{
  xfree(b->bitmap);
  xfree(b->patterns);
  xfree(b);
}
