/*
    SWARM

    Copyright (C) 2012-2025 Torbjorn Rognes and Frederic Mahe

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
#include <vector>


auto hash_getindex(uint64_t hash) -> uint64_t;

auto hash_getnextindex(uint64_t index) -> uint64_t;

auto hash_set_occupied(uint64_t index) -> void;

auto hash_is_occupied(uint64_t index) -> bool;

auto hash_set_value(uint64_t index, uint64_t hash) -> void;

auto hash_compare_value(uint64_t index, uint64_t hash) -> bool;

auto hash_get_data(uint64_t index) -> unsigned int;

auto hash_set_data(uint64_t index, unsigned int amplicon_id) -> void;

auto hash_alloc(uint64_t amplicons,
                std::vector<unsigned char>& hash_occupied_v,
                std::vector<uint64_t>& hash_values_v,
                std::vector<unsigned int>& hash_data_v) -> uint64_t;

auto hash_free() -> void;
