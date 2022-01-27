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

#include <cstdint>

extern uint64_t hash_mask;
extern unsigned char * hash_occupied;
extern uint64_t * hash_values;
extern unsigned int * hash_data;


inline auto hash_getindex(uint64_t hash) -> uint64_t
{
  // Shift bits right to get independence from the simple Bloom filter hash
  constexpr unsigned int divider {32};  // drop the first 32 bits
  hash = hash >> divider;
  return hash & hash_mask;
}

inline auto hash_getnextindex(uint64_t j) -> uint64_t
{
  return (j+1) & hash_mask;
}

inline void hash_set_occupied(uint64_t j)
{
  constexpr unsigned int divider {3};  // drop the first 3 bits
  constexpr unsigned int max_range {7};  // j & max_range = values ranging from 0 to 7
  hash_occupied[j >> divider] |= (1 << (j & max_range));
}

inline auto hash_is_occupied(uint64_t j) -> bool
{
  constexpr unsigned int divider {3};  // drop the first 3 bits
  constexpr unsigned int max_range {7};  // j & max_range = values ranging from 0 to 7
  return (hash_occupied[j >> divider] & (1 << (j & max_range))) != 0;
}

inline void hash_set_value(uint64_t j, uint64_t hash)
{
  hash_values[j] = hash;
}

inline auto hash_compare_value(uint64_t j, uint64_t hash) -> bool
{
  return (hash_values[j] == hash);
}

inline auto hash_get_data(uint64_t j) -> unsigned int
{
  return hash_data[j];
}

inline void hash_set_data(uint64_t j, unsigned int x)
{
  hash_data[j] = x;
}

void hash_zap(uint64_t hashtablesize);

auto hash_alloc(uint64_t amplicons) -> uint64_t;

void hash_free();
