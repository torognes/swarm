/*
    SWARM

    Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

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

#include <algorithm>  // std::fill
#include <cassert>
#include <cstdint>
#include "hashtable.h"
#include "hashtable_size.h"


auto Hashtable::allocate(const uint64_t amplicons) -> uint64_t
{
  static constexpr int padding {63};  // make sure our final value is >= 64 / 8
  static constexpr int convert_to_bytes {8};

  const auto hashtablesize = compute_hashtable_size(amplicons);
  mask = hashtablesize - 1;

  occupied.assign((hashtablesize + padding) / convert_to_bytes, 0U);
  values.assign(hashtablesize, 0U);
  data.assign(hashtablesize, 0U);

  return hashtablesize;
}


auto Hashtable::clear() -> void
{
  std::fill(occupied.begin(), occupied.end(), 0U);
}


auto Hashtable::getindex(uint64_t hash) const noexcept -> uint64_t
{
  // Shift bits right to get independence from the simple Bloom filter hash
  static constexpr auto divider = 32U;  // drop the first 32 bits
  hash = hash >> divider;
  return hash & mask;
}


auto Hashtable::getnextindex(const uint64_t index) const noexcept -> uint64_t
{
  return (index + 1) & mask;
}


auto Hashtable::set_occupied(const uint64_t index) noexcept -> void
{
  static constexpr auto divider = 3U;
  static constexpr auto max_range = 7U;  // 0000 0111
  auto const multiplier = index & max_range;  // mask all but the first 3 bits
  assert(multiplier <= 7);
  auto const bit_to_set = static_cast<unsigned char>(1U << multiplier);  // bit 0 to bit 7
  auto const position = index >> divider;  // divide by 8, so drop the first 3 bits
  assert(position < occupied.size());
  occupied[position] |= bit_to_set;
}


auto Hashtable::is_occupied(const uint64_t index) const noexcept -> bool
{
  static constexpr auto divider = 3U;
  static constexpr auto max_range = 7U;
  auto const multiplier = index & max_range;  // mask all but the first 3 bits
  assert(multiplier <= 7);
  auto const bit_to_check = static_cast<unsigned char>(1U << multiplier);  // bit 0 to bit 7
  auto const position = index >> divider;  // divide by 8, so drop the first 3 bits
  assert(position < occupied.size());
  return (occupied[position] & bit_to_check) != 0;
}


auto Hashtable::set_value(const uint64_t index, const uint64_t hash) noexcept -> void
{
  assert(index < values.size());
  values[index] = hash;
}


auto Hashtable::compare_value(const uint64_t index, const uint64_t hash) const noexcept -> bool
{
  assert(index < values.size());
  return values[index] == hash;
}


auto Hashtable::get_data(const uint64_t index) const noexcept -> unsigned int
{
  assert(index < data.size());
  return data[index];
}


auto Hashtable::set_data(const uint64_t index, const unsigned int amplicon_id) noexcept -> void
{
  assert(index < data.size());
  data[index] = amplicon_id;
}
