/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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
#include <cstdio>  // size_t


auto xmalloc(std::size_t size) -> void *;
auto xrealloc(void * ptr, std::size_t size) -> void *;
auto xfree(void * ptr) -> void;
auto xgetline(char ** linep, std::size_t * linecapp, FILE * stream) -> ssize_t;
auto progress_init(const char * prompt, uint64_t size) -> void;
auto progress_update(uint64_t progress) -> void;
auto progress_done() -> void;
auto fopen_input(const char * filename) -> std::FILE *;
auto fopen_output(const char * filename) -> std::FILE *;
