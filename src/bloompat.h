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

constexpr unsigned int BLOOM_PATTERN_SHIFT {10};
constexpr unsigned int BLOOM_PATTERN_COUNT {1 << BLOOM_PATTERN_SHIFT};
constexpr unsigned int BLOOM_PATTERN_MASK {BLOOM_PATTERN_COUNT - 1};

struct bloom_s
{
  uint64_t size;
  uint64_t mask;
  uint64_t * bitmap;
  uint64_t patterns[BLOOM_PATTERN_COUNT];
};

void bloom_zap(struct bloom_s * b);

struct bloom_s * bloom_init(uint64_t size);

void bloom_exit(struct bloom_s * b);

inline uint64_t * bloom_adr(struct bloom_s * b, uint64_t h)
{
  return b->bitmap + ((h >> BLOOM_PATTERN_SHIFT) & b->mask);
}

inline uint64_t bloom_pat(struct bloom_s * b, uint64_t h)
{
  return b->patterns[h & BLOOM_PATTERN_MASK];
}

inline void bloom_set(struct bloom_s * b, uint64_t h)
{
  * bloom_adr(b, h) &= ~ bloom_pat(b, h);
}

inline bool bloom_get(struct bloom_s * b, uint64_t h)
{
  return ! (* bloom_adr(b, h) & bloom_pat(b, h));
}
