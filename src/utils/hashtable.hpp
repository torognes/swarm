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

#ifndef SWARM_UTILS_HASHTABLE_H
#define SWARM_UTILS_HASHTABLE_H

#include <cstdint>
#include <vector>


// Open-addressing hash table with linear probing, specialised for the
// (64-bit Zobrist hash -> amplicon id) lookups performed by algod1.cpp.
// Three flat buffers back the table: a packed bitset of occupancy, the
// stored hash value at each slot, and the amplicon id at each slot.
class Hashtable {
public:

  // Allocate buffers sized for `amplicons` entries; returns the resulting
  // hash-table size (always a power of two).
  auto allocate(uint64_t amplicons) -> uint64_t;

  // Reset all occupancy bits to zero, preserving allocated capacity. The
  // value/data buffers are left untouched: they are only read when the
  // matching occupancy bit is set.
  auto clear() -> void;

  auto getindex(uint64_t hash) const noexcept -> uint64_t;
  auto getnextindex(uint64_t index) const noexcept -> uint64_t;

  auto set_occupied(uint64_t index) noexcept -> void;
  auto is_occupied(uint64_t index) const noexcept -> bool;

  auto set_value(uint64_t index, uint64_t hash) noexcept -> void;
  auto compare_value(uint64_t index, uint64_t hash) const noexcept -> bool;

  auto get_data(uint64_t index) const noexcept -> unsigned int;
  auto set_data(uint64_t index, unsigned int amplicon_id) noexcept -> void;

private:

  uint64_t mask {0};
  std::vector<unsigned char> occupied;
  std::vector<uint64_t> values;
  std::vector<unsigned int> data;
};

#endif  // SWARM_UTILS_HASHTABLE_H
