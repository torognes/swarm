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

struct bloomflex_s
{
  uint64_t size; /* size in number of longs (8 bytes) */
  uint64_t pattern_shift;
  uint64_t pattern_count;
  uint64_t pattern_mask;
  uint64_t pattern_k;
  uint64_t * bitmap;
  uint64_t * patterns;
};

struct bloomflex_s * bloomflex_init(uint64_t size, unsigned int k);

void bloomflex_exit(struct bloomflex_s * b);

inline uint64_t * bloomflex_adr(struct bloomflex_s * b, uint64_t h)
{
  return b->bitmap + ((h >> b->pattern_shift) % b->size);
}

inline uint64_t bloomflex_pat(struct bloomflex_s * b,
                                     uint64_t h)
{
  return b->patterns[h & b->pattern_mask];
}

inline void bloomflex_set(struct bloomflex_s * b, uint64_t h)
{
  * bloomflex_adr(b, h) &= ~ bloomflex_pat(b, h);
}

inline bool bloomflex_get(struct bloomflex_s * b, uint64_t h)
{
  return ! (* bloomflex_adr(b, h) & bloomflex_pat(b, h));
}
