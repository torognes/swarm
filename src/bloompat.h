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

#include <array>
#include <cstdint>  // uint64_t
#include <vector>

constexpr unsigned int bloom_pattern_shift {10};
constexpr unsigned int bloom_pattern_count {1U << bloom_pattern_shift};
constexpr unsigned int bloom_pattern_mask {bloom_pattern_count - 1};

struct bloom_s
{
  uint64_t size;
  uint64_t mask;
  std::vector<uint64_t> bitmap_v;
  uint64_t * bitmap;
  std::array<uint64_t, bloom_pattern_count> patterns;
};

auto bloom_zap(struct bloom_s & bloom_filter) -> void;

auto bloom_init(uint64_t size, struct bloom_s & bloom_filter) -> struct bloom_s *;

auto bloom_exit(struct bloom_s * bloom_filter) -> void;


// --------------------------------------------------- used only in this header

inline auto bloom_adr(struct bloom_s * bloom_filter, uint64_t h) -> uint64_t *
{
  return bloom_filter->bitmap + ((h >> bloom_pattern_shift) & bloom_filter->mask);
}

inline auto bloom_pat(struct bloom_s * bloom_filter, uint64_t h) -> uint64_t
{
  return bloom_filter->patterns[h & bloom_pattern_mask];
}

// ------------------------------------------------------------------------ end

// used in algod1.cc
inline auto bloom_set(struct bloom_s * bloom_filter, uint64_t h) -> void
{
  * bloom_adr(bloom_filter, h) &= compl bloom_pat(bloom_filter, h);
}

// used in algod1.cc
inline auto bloom_get(struct bloom_s * bloom_filter, uint64_t h) -> bool
{
  return (* bloom_adr(bloom_filter, h) & bloom_pat(bloom_filter, h)) == 0U;
}
