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

#include <cstdint>  // uint64_t
#include <vector>


struct bloomflex_s
{
  uint64_t size = 0; /* size in number of longs (8 bytes) */
  uint64_t pattern_shift = 0;
  uint64_t pattern_count = 0;
  uint64_t pattern_mask = 0;
  uint64_t pattern_k = 0;
  std::vector<uint64_t> bitmap_v;
  uint64_t * bitmap = nullptr;
  std::vector<uint64_t> patterns_v;
  uint64_t * patterns = nullptr;
};

auto bloomflex_init(uint64_t size, unsigned int n_hash_functions,
                    struct bloomflex_s & bloom_filter) -> struct bloomflex_s *;

auto bloomflex_exit(struct bloomflex_s & bloom_filter) -> void;

auto bloomflex_adr(struct bloomflex_s * bloom_filter, uint64_t hash) -> uint64_t *;

auto bloomflex_pat(struct bloomflex_s * bloom_filter, uint64_t hash) -> uint64_t;

auto bloomflex_set(struct bloomflex_s * bloom_filter, uint64_t hash) -> void;

auto bloomflex_get(struct bloomflex_s * bloom_filter, uint64_t hash) -> bool;
