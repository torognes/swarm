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

uint64_t * zobrist_tab_base = 0;

void zobrist_init(unsigned int n)
{
  /*
    Generate n random 64-bit numbers, where n is 4 * longest.
    They will represent the four different bases in any position (0-1023)
    in the sequence. They will be XOR'ed together to form the hash of
    that sequence.
  */

  zobrist_tab_base = (uint64_t *) xmalloc(8 * 4 * n);

  for (unsigned int i = 0; i < 4 * n; i++)
    {
      uint64_t z;
      z = arch_random();
      z <<= 16;
      z ^= arch_random();
      z <<= 16;
      z ^= arch_random();
      z <<= 16;
      z ^= arch_random();
      zobrist_tab_base[i] = z;
    }
}

void zobrist_exit()
{
  xfree(zobrist_tab_base);
}

uint64_t zobrist_hash(unsigned char * s, unsigned int len)
{
  /* compute the Zobrist hash function of sequence s of length len. */
  /* len is the actual number of bases in the sequence */
  /* it is encoded in (len+3)/4 bytes */

  uint64_t * q = (uint64_t *) s;
  uint64_t x = 0;
  uint64_t z = 0;
  for(unsigned int p = 0; p < len; p++)
    {
      if ((p & 31) == 0)
        x = q[p / 32];
      else
        x >>= 2;
      z ^= zobrist_value(p, x & 3);
    }
  return z;
}

uint64_t zobrist_hash_delete_first(unsigned char * s, unsigned int len)
{
  /* compute the Zobrist hash function of sequence s,
     but delete the first base */

  uint64_t * q = (uint64_t *) s;
  uint64_t x = q[0];
  uint64_t z = 0;
  for(unsigned int p = 1; p < len; p++)
    {
      if ((p & 31) == 0)
        x = q[p / 32];
      else
        x >>= 2;
      z ^= zobrist_value(p - 1, x & 3);
    }
  return z;
}

uint64_t zobrist_hash_insert_first(unsigned char * s, unsigned int len)
{
  /* compute the Zobrist hash function of sequence s,
     but insert a gap (no value) before the first base */

  uint64_t * q = (uint64_t *) s;
  uint64_t x = 0;
  uint64_t z = 0;
  for(unsigned int p = 0; p < len; p++)
    {
      if ((p & 31) == 0)
        x = q[p / 32];
      else
        x >>= 2;
      z ^= zobrist_value(p + 1, x & 3);
    }
  return z;
}
