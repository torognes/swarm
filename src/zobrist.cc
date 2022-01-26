/*
    SWARM

    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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
#include "pseudo_rng.h"
#include "util.h"
#include "zobrist.h"


uint64_t * zobrist_tab_base = nullptr;
uint64_t * zobrist_tab_byte_base = nullptr;


void zobrist_init(const unsigned int n)
{
  /*
    Generate 4n random 64-bit numbers. They will represent the four
    different bases in any position (1 to n) of a sequence.
    They will be XOR'ed together to form the hash of that sequence.
    The number n should be the length of the longest sequence to be
    hashed including potential additional insertions.

    The number is generated by xor'ing together four shifted
    31-bit random numbers.
  */
  constexpr unsigned int byte_range {256};
  constexpr unsigned int multiplier {16};

  /* allocate memory for tables */

  zobrist_tab_base = static_cast<uint64_t *>
    (xmalloc(4 * n * sizeof(uint64_t)));

  zobrist_tab_byte_base = static_cast<uint64_t *>
    (xmalloc(byte_range * (n / 4) * sizeof(uint64_t)));

  /* fill table with random 64 bit numbers */

  for(auto i = 0U; i < 4 * n; i++)
    {
      auto z = 0ULL;
      z = rand_64();
      z <<= multiplier;
      z ^= rand_64();
      z <<= multiplier;
      z ^= rand_64();
      z <<= multiplier;
      z ^= rand_64();
      zobrist_tab_base[i] = z;
    }

  /* combine into bytes for faster computations */

  for(auto i = 0U; i < n / 4; i++) {
    for(auto j = 0U; j < byte_range; j++) {
      auto z = 0ULL;
      auto x = j;
      z ^= zobrist_value(4 * i + 0, x & 3);
      x >>= 2;
      z ^= zobrist_value(4 * i + 1, x & 3);
      x >>= 2;
      z ^= zobrist_value(4 * i + 2, x & 3);
      x >>= 2;
      z ^= zobrist_value(4 * i + 3, x & 3);
      zobrist_tab_byte_base[byte_range * i + j] = z;
    }
  }
}

void zobrist_exit()
{
  xfree(zobrist_tab_byte_base);
  xfree(zobrist_tab_base);
}


auto zobrist_hash(unsigned char * seq, unsigned int len) -> uint64_t
{
  /* compute the Zobrist hash function of sequence seq of length len. */
  /* len is the actual number of bases in the sequence */
  /* it is encoded in (len + 3 ) / 4 bytes */

  constexpr unsigned int offset {64};
  constexpr unsigned int nt_per_uint64 {32};  // 32 nucleotides can fit in a uint64
  auto * q = reinterpret_cast<uint64_t *>(seq);
  uint64_t z = 0;
  unsigned int p = 0;
  auto * qb = reinterpret_cast<unsigned char *>(q);

  while (p + nt_per_uint64 < len)
    {
      for(auto i = 0U; i < nt_per_uint64; i += 4) {
        // i = {0, 4, 8, 12, 16, 20, 24, 28}
        z ^= zobrist_tab_byte_base[offset * (p + i) + *qb++];
      }
      p += nt_per_uint64;
    }

  while (p + 4 < len)
    {
      z ^= zobrist_tab_byte_base[offset * p + *qb++];
      p += 4;
    }

  if (p < len)
    {
      uint64_t x = *qb++;
      while (p < len)
        {
          z ^= zobrist_value(p, x & 3);
          x >>= 2;
          p++;
        }
    }

  return z;
}


auto zobrist_hash_delete_first(unsigned char * seq, unsigned int len) -> uint64_t
{
  /* compute the Zobrist hash function of sequence seq,
     but delete the first base */

  constexpr unsigned int nt_per_uint64 {32};  // 32 nucleotides can fit in a uint64
  auto * q = reinterpret_cast<uint64_t *>(seq);
  uint64_t x = q[0];
  uint64_t z = 0;
  for(auto p = 1U; p < len; p++)
    {
      if ((p & (nt_per_uint64 - 1)) == 0) {
        x = q[p / nt_per_uint64];
      }
      else {
        x >>= 2;
      }
      z ^= zobrist_value(p - 1, x & 3);
    }
  return z;
}

auto zobrist_hash_insert_first(unsigned char * seq, unsigned int len) -> uint64_t
{
  /* compute the Zobrist hash function of sequence seq,
     but insert a gap (no value) before the first base */

  constexpr unsigned int nt_per_uint64 {32};  // 32 nucleotides can fit in a uint64
  auto * q = reinterpret_cast<uint64_t *>(seq);
  uint64_t x = 0;
  uint64_t z = 0;
  for(auto p = 0U; p < len; p++)
    {
      if ((p & (nt_per_uint64 - 1)) == 0) {
        x = q[p / nt_per_uint64];
      }
      else {
        x >>= 2;
      }
      z ^= zobrist_value(p + 1, x & 3);
    }
  return z;
}
