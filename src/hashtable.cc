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
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>
#include <cstring>  // memset
#include <iterator>  // std::next
#include <vector>
#include "hashtable.h"
#include "utils/hashtable_size.h"


uint64_t hash_mask {0};
unsigned char * hash_occupied {nullptr};
uint64_t * hash_values {nullptr};
unsigned int * hash_data {nullptr};


auto hash_getindex(uint64_t hash) -> uint64_t
{
  // Shift bits right to get independence from the simple Bloom filter hash
  static constexpr auto divider = 32U;  // drop the first 32 bits
  hash = hash >> divider;
  return hash & hash_mask;
}


auto hash_getnextindex(uint64_t index) -> uint64_t
{
  return (index + 1) & hash_mask;
}


auto hash_set_occupied(const uint64_t index) -> void
{
  static constexpr auto divider = 3U;
  static constexpr auto max_range = 7U;
  auto const multiplier = index & max_range;  // mask all but the first 3 bits
  assert(multiplier <= 7);
  auto const bit_to_set = static_cast<unsigned char>(1U << multiplier);  // bit 0 to bit 7
  auto const position = static_cast<std::ptrdiff_t>(index >> divider);  // divide by 8, so drop the first 3 bits
  auto & target_bucket = *std::next(hash_occupied, position);
  target_bucket |= bit_to_set;
}


auto hash_is_occupied(const uint64_t index) -> bool
{
  static constexpr auto divider = 3U;
  static constexpr auto max_range = 7U;
  auto const multiplier = index & max_range;  // mask all but the first 3 bits
  assert(multiplier <= 7);
  auto const bit_to_check = static_cast<unsigned char>(1U << multiplier);  // bit 0 to bit 7
  auto const position = static_cast<std::ptrdiff_t>(index >> divider);  // divide by 8, so drop the first 3 bits
  auto & target_bucket = *std::next(hash_occupied, position);
  return (target_bucket & bit_to_check) != 0;
}


auto hash_set_value(uint64_t index, uint64_t hash) -> void
{
  hash_values[index] = hash;
}


auto hash_compare_value(uint64_t index, uint64_t hash) -> bool
{
  return (hash_values[index] == hash);
}


auto hash_get_data(uint64_t index) -> unsigned int
{
  return hash_data[index];
}


auto hash_set_data(uint64_t index, unsigned int amplicon_id) -> void
{
  hash_data[index] = amplicon_id;
}


auto hash_alloc(const uint64_t amplicons,
                std::vector<unsigned char>& hash_occupied_v,
                std::vector<uint64_t>& hash_values_v,
                std::vector<unsigned int>& hash_data_v) -> uint64_t
{
  static constexpr int padding {63};  // make sure our final value is >= 64 / 8
  static constexpr int convert_to_bytes {8};

  const auto hashtablesize {compute_hashtable_size(amplicons)};
  hash_mask = hashtablesize - 1;

  hash_occupied_v.resize((hashtablesize + padding) / convert_to_bytes);
  hash_occupied = hash_occupied_v.data();

  hash_values_v.resize(hashtablesize);
  hash_values = hash_values_v.data();

  hash_data_v.resize(hashtablesize);
  hash_data = hash_data_v.data();

  return hashtablesize;
}


auto hash_free() -> void
{
  hash_occupied = nullptr;
  hash_values = nullptr;
  hash_data = nullptr;
}
