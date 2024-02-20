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


// refactoring: use a struct to communicate and avoid global variables
// struct Progress_status {
//   const char * progress_prompt = nullptr;
//   uint64_t progress_next = 0;
//   uint64_t progress_size = 0;
//   uint64_t progress_chunk = 0;
// };

auto progress_init(const char * prompt, uint64_t size) -> void;
auto progress_update(uint64_t progress) -> void;
auto progress_done(struct Parameters const & parameters) -> void;
