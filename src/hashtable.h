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

#include <cassert>
#include <cstdint>
#include <vector>


extern uint64_t hash_mask;
extern unsigned char * hash_occupied;
extern uint64_t * hash_values;
extern unsigned int * hash_data;


inline auto hash_getindex(uint64_t hash) -> uint64_t
{
  // Shift bits right to get independence from the simple Bloom filter hash
  static constexpr auto divider = 32U;  // drop the first 32 bits
  hash = hash >> divider;
  return hash & hash_mask;
}

inline auto hash_getnextindex(uint64_t index) -> uint64_t
{
  return (index + 1) & hash_mask;
}

inline auto hash_set_occupied(uint64_t index) -> void
{
  static constexpr auto divider = 3U;  // drop the first 3 bits
  static constexpr auto max_range = 7U;  // index & max_range = values ranging from 0 to 7
  assert((index & max_range) <= 7);
  hash_occupied[index >> divider] |= static_cast<unsigned char>((1U << (index & max_range)));
}

inline auto hash_is_occupied(uint64_t index) -> bool
{
  static constexpr auto divider = 3U;  // drop the first 3 bits
  static constexpr auto max_range = 7U;  // index & max_range = values ranging from 0 to 7
  return (hash_occupied[index >> divider] & (1U << (index & max_range))) != 0;
}

inline auto hash_set_value(uint64_t index, uint64_t hash) -> void
{
  hash_values[index] = hash;
}

inline auto hash_compare_value(uint64_t index, uint64_t hash) -> bool
{
  return (hash_values[index] == hash);
}

inline auto hash_get_data(uint64_t index) -> unsigned int
{
  return hash_data[index];
}

inline auto hash_set_data(uint64_t index, unsigned int amplicon_id) -> void
{
  hash_data[index] = amplicon_id;
}

auto hash_alloc(uint64_t amplicons,
                std::vector<unsigned char>& hash_occupied_v,
                std::vector<uint64_t>& hash_values_v,
                std::vector<unsigned int>& hash_data_v) -> uint64_t;

auto hash_free() -> void;
