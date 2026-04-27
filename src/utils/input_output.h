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

#ifndef SWARM_UTILS_INPUT_OUTPUT_H
#define SWARM_UTILS_INPUT_OUTPUT_H

#include <cstdio>  // FILE, fclose
#include <memory>  // unique_ptr


// RAII wrapper for std::FILE *: the deleter calls std::fclose, so a
// FileHandle that goes out of scope automatically closes its file.
//
// Note: taking the address of a standard library function (such as
// &std::fclose) as deleter is unspecified behaviour; prefer a deleter
// struct with an operator() that calls std::fclose.
struct CloseFileHandle {
  auto operator()(std::FILE * file_handle) const -> void {
    static_cast<void>(std::fclose(file_handle));
  }
};

using FileHandle = std::unique_ptr<std::FILE, CloseFileHandle>;


auto fopen_input(const char * filename) -> FileHandle;
auto fopen_output(const char * filename) -> FileHandle;

#endif  // SWARM_UTILS_INPUT_OUTPUT_H
