/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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

#include <algorithm>  // std::for_each
#include <cstdint>  // uint64_t
#include <vector>
#include "utils/pseudo_rng.h"
#include "zobrist.h"


uint64_t * zobrist_tab_base = nullptr;
uint64_t * zobrist_tab_byte_base = nullptr;


auto fill_rng_table(const unsigned int zobrist_len,
                    std::vector<uint64_t> & zobrist_tab_base_v) -> void
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
  static constexpr auto multiplier = 16U;

  /* allocate base table and fill with random 64 bit numbers */
  zobrist_tab_base_v.resize(4 * zobrist_len);
  zobrist_tab_base = zobrist_tab_base_v.data();

  std::for_each(zobrist_tab_base_v.begin(),
                zobrist_tab_base_v.end(),
                [](uint64_t & rng_value) {
                  // refactoring: comment states 31-bit random numbers?!
                  rng_value = rand_64();
                  rng_value <<= multiplier;
                  rng_value ^= rand_64();
                  rng_value <<= multiplier;
                  rng_value ^= rand_64();
                  rng_value <<= multiplier;
                  rng_value ^= rand_64();
                });
}


auto fill_rng_byte_table(const unsigned int zobrist_len,
                         std::vector<uint64_t> const & zobrist_tab_base_v,
                         std::vector<uint64_t> & zobrist_tab_byte_base_v) -> void
{
  static constexpr auto byte_range = 256U;

  /* allocate byte table and combine into bytes for faster computations */
  zobrist_tab_byte_base_v.resize(byte_range * (zobrist_len / 4));
  zobrist_tab_byte_base = zobrist_tab_byte_base_v.data();

  for(auto i = 0U; i < zobrist_len / 4; i++) {
    for(auto j = 0U; j < byte_range; j++) {
      auto rng_value = 0ULL;
      auto offset = j;
      // rng value stored at: 4 *  position   +  offset & 3U (= 0, 1, 2, or 3)
      rng_value ^= zobrist_tab_base_v[4 * (4 * i + 0) + (offset & 3U)];
      offset >>= 2U;
      rng_value ^= zobrist_tab_base_v[4 * (4 * i + 1) + (offset & 3U)];
      offset >>= 2U;
      rng_value ^= zobrist_tab_base_v[4 * (4 * i + 2) + (offset & 3U)];
      offset >>= 2U;
      rng_value ^= zobrist_tab_base_v[4 * (4 * i + 3) + (offset & 3U)];
      zobrist_tab_byte_base_v[byte_range * i + j] = rng_value;
    }
  }
}


auto zobrist_init(const unsigned int zobrist_len,
                  std::vector<uint64_t> & zobrist_tab_base_v,
                  std::vector<uint64_t> & zobrist_tab_byte_base_v) -> void
{
  fill_rng_table(zobrist_len, zobrist_tab_base_v);
  fill_rng_byte_table(zobrist_len, zobrist_tab_base_v, zobrist_tab_byte_base_v);
}


auto zobrist_exit() -> void
{
  zobrist_tab_byte_base = nullptr;
  zobrist_tab_base = nullptr;
}


auto zobrist_hash(unsigned char * seq, const unsigned int len) -> uint64_t
{
  /* compute the Zobrist hash function of sequence seq of length len. */
  /* len is the actual number of bases in the sequence */
  /* it is encoded in (len + 3 ) / 4 bytes */

  static constexpr auto offset = 64U;
  static constexpr auto nt_per_uint64 = 32U;  // 32 nucleotides can fit in a uint64
  // refactoring: why interpret first as uint64_t then as uint8_t??
  auto * query = reinterpret_cast<uint64_t *>(seq);
  uint64_t zobrist_hash = 0;
  auto pos = 0U;
  auto * query_in_bytes = reinterpret_cast<unsigned char *>(query);

  while (pos + nt_per_uint64 < len)
    {
      for(auto i = 0U; i < nt_per_uint64; i += 4) {
        // i = {0, 4, 8, 12, 16, 20, 24, 28}
        zobrist_hash ^= zobrist_tab_byte_base[offset * (pos + i) + *query_in_bytes];
        ++query_in_bytes;
      }
      pos += nt_per_uint64;
    }

  while (pos + 4 < len)
    {
      zobrist_hash ^= zobrist_tab_byte_base[offset * pos + *query_in_bytes];
      ++query_in_bytes;
      pos += 4;
    }

  if (pos < len)
    {
      uint64_t next_byte = *query_in_bytes;
      ++query_in_bytes;
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

  static constexpr auto nt_per_uint64 = 32U;  // 32 nucleotides can fit in a uint64
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
