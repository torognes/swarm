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

#include "swarm.h"

unsigned long * zobrist_tab_base = 0;

static unsigned long * zobrist_tab_byte = 0;

void zobrist_init(unsigned int longest)
{
  /*
    Generate n random 64-bit numbers, where n is 4 * longest.
    They will represent the four different bases in any position (0-1023)
    in the sequence. They will be XOR'ed together to form the hash of
    that sequence.
  */

  srandom(1);

  zobrist_tab_base = (unsigned long *) xmalloc(8 * 4 * longest);

  for (unsigned int i = 0; i < 4 * longest; i++)
    {
      unsigned long z;
      z = random();
      z <<= 16;
      z ^= random();
      z <<= 16;
      z ^= random();
      z <<= 16;
      z ^= random();
      zobrist_tab_base[i] = z;
    }


  /* 
     For increased performance, we make a table of precomputed values
     for every possible 4-base sequence (represented as one byte) in any
     position (0, 4, 8, ...).
  */

  zobrist_tab_byte = (unsigned long *) xmalloc(8 * 256 * (longest/4));

  for (unsigned int p = 0; p < longest / 4; p++)
    for (unsigned int x = 0; x < 256; x++)
      {
        unsigned long z = 0;
        for (unsigned int i = 0; i < 4; i++)
          z ^= zobrist_tab_base[(4 * (4 * p + i)) + ((x >> (2*i)) & 3)];
        zobrist_tab_byte[256 * p + x] = z;

#if 0
        fprintf(stderr,
                "Zobrist hash pos %3u byte %02x: %016lx\n",
                p,
                x,
                z);
#endif
      }
}

void zobrist_exit()
{
  free(zobrist_tab_byte);
  free(zobrist_tab_base);
}

unsigned long zobrist_hash(unsigned char * s, unsigned int len)
{
  /* compute the Zobrist hash function of sequence s of length len. */
  /* len is the actual number of bases in the sequence */
  /* it is encoded in (len+3)/4 bytes */

  unsigned long z = 0;

  for (unsigned int i = 0; i < (len / 4); i++)
    z ^= zobrist_tab_byte[256 * i + (unsigned int)(s[i])];

  if ((len & 3) > 0)
    {
      unsigned char x = s[len/4];
      for(unsigned int p = 4 * (len / 4); p < len; p++)
        {
          z ^= zobrist_value(p, x & 3);
          x >>= 2;
        }
    }

#if 0
  unsigned long z2 = 0;
  for (unsigned int p = 0; p < len ; p++)
    {
      unsigned char x = (s[p/4] >> (2*(p & 3))) & 3;
      z2 ^= zobrist_tab_base[4*p + x];
    }
  if (z != z2)
    {
      fprintf(stderr,
              "Zobrist mismatch: z=%016lx z2=%016lx\n",
              z,
              z2);
    }
#endif

  return z;
}
