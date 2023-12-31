/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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

#include <cstdint>
#include "utils/pseudo_rng.h"
#include "zobrist.h"


uint64_t * zobrist_tab_base = nullptr;
uint64_t * zobrist_tab_byte_base = nullptr;


void zobrist_init(const unsigned int zobrist_len)
{
  /*
    Generate 4n random 64-bit numbers. They will represent the four
    different bases in any position (1 to n) of a sequence.  They will
    be XOR'ed together to form the hash of that sequence.  The number
    n, here 'zobrist_len', should be the length of the longest
    sequence to be hashed including potential additional insertions.

    The number is generated by xor'ing together four shifted
    31-bit random numbers.
  */
  static constexpr unsigned int byte_range {256};
  static constexpr unsigned int multiplier {16};

  /* allocate base table and fill with random 64 bit numbers */

  zobrist_tab_base = new uint64_t[4 * zobrist_len];

  for(auto i = 0U; i < 4 * zobrist_len; i++)
    {
      auto z = 0ULL;
      z = rand_64();  // refactoring: comment states 31-bit random numbers?!
      z <<= multiplier;
      z ^= rand_64();
      z <<= multiplier;
      z ^= rand_64();
      z <<= multiplier;
      z ^= rand_64();
      zobrist_tab_base[i] = z;
    }

  /* allocate byte table and combine into bytes for faster computations */

  zobrist_tab_byte_base = new uint64_t[byte_range * (zobrist_len / 4)];

  for(auto i = 0U; i < zobrist_len / 4; i++) {
    for(auto j = 0U; j < byte_range; j++) {
      auto z = 0ULL;
      auto x = j;
      z ^= zobrist_value(4 * i + 0, x & 3U);
      x >>= 2U;
      z ^= zobrist_value(4 * i + 1, x & 3U);
      x >>= 2U;
      z ^= zobrist_value(4 * i + 2, x & 3U);
      x >>= 2U;
      z ^= zobrist_value(4 * i + 3, x & 3U);
      zobrist_tab_byte_base[byte_range * i + j] = z;
    }
  }
}

void zobrist_exit()
{
  delete [] zobrist_tab_byte_base;
  zobrist_tab_byte_base = nullptr;
  delete [] zobrist_tab_base;
  zobrist_tab_base = nullptr;
}


auto zobrist_hash(unsigned char * seq, const unsigned int len) -> uint64_t
{
  /* compute the Zobrist hash function of sequence seq of length len. */
  /* len is the actual number of bases in the sequence */
  /* it is encoded in (len + 3 ) / 4 bytes */

  static constexpr unsigned int offset {64};
  static constexpr unsigned int nt_per_uint64 {32};  // 32 nucleotides can fit in a uint64
  auto * query = reinterpret_cast<uint64_t *>(seq);
  uint64_t zobrist_hash = 0;
  unsigned int pos = 0;
  auto * query_in_bytes = reinterpret_cast<unsigned char *>(query);

  while (pos + nt_per_uint64 < len)
    {
      for(auto i = 0U; i < nt_per_uint64; i += 4) {
        // i = {0, 4, 8, 12, 16, 20, 24, 28}
        zobrist_hash ^= zobrist_tab_byte_base[offset * (pos + i) + *query_in_bytes++];
      }
      pos += nt_per_uint64;
    }

  while (pos + 4 < len)
    {
      zobrist_hash ^= zobrist_tab_byte_base[offset * pos + *query_in_bytes++];
      pos += 4;
    }

  if (pos < len)
    {
      uint64_t next_byte = *query_in_bytes++;
      while (pos < len)
        {
          zobrist_hash ^= zobrist_value(pos, next_byte & 3U);
          next_byte >>= 2U;
          ++pos;
        }
    }

  return zobrist_hash;
}


auto zobrist_hash_delete_first(unsigned char * seq, const unsigned int len) -> uint64_t
{
  /* compute the Zobrist hash function of sequence seq,
     but delete the first base */

  static constexpr unsigned int nt_per_uint64 {32};  // 32 nucleotides can fit in a uint64
  auto * q = reinterpret_cast<uint64_t *>(seq);
  uint64_t x = q[0];
  uint64_t zobrist_hash = 0;
  for(auto pos = 1U; pos < len; pos++)
    {
      if ((pos & (nt_per_uint64 - 1)) == 0) {
        x = q[pos / nt_per_uint64];  // UBSAN: misaligned address for type 'long unsigned int', which requires 8 byte alignment
      }
      else {
        x >>= 2U;
      }
      zobrist_hash ^= zobrist_value(pos - 1, x & 3U);
    }
  return zobrist_hash;
}

auto zobrist_hash_insert_first(unsigned char * seq, const unsigned int len) -> uint64_t
{
  /* compute the Zobrist hash function of sequence seq,
     but insert a gap (no value) before the first base */

  static constexpr unsigned int nt_per_uint64 {32};  // 32 nucleotides can fit in a uint64
  auto * q = reinterpret_cast<uint64_t *>(seq);
  uint64_t x = 0;
  uint64_t zobrist_hash = 0;
  for(auto pos = 0U; pos < len; pos++)
    {
      if ((pos & (nt_per_uint64 - 1)) == 0) {
        x = q[pos / nt_per_uint64];  // UBSAN: misaligned address for type 'long unsigned int', which requires 8 byte alignment
      }
      else {
        x >>= 2U;
      }
      zobrist_hash ^= zobrist_value(pos + 1, x & 3U);
    }
  return zobrist_hash;
}
