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

#include <cstdio>  // FILE
#include <memory>  // unique_ptr
#include <string>  // std::string


// RAII wrapper for std::FILE *: the deleter calls std::fclose, so a
// FileHandle that goes out of scope automatically closes its file.
//
// Note: taking the address of a standard library function (such as
// &std::fclose) as deleter is unspecified behaviour; prefer a deleter
// struct with an operator() that calls std::fclose.
// The close is also where I/O failures are noticed. Nothing else in swarm
// checks a write: every fprint/fputc/fwrite return value is discarded, so
// before this a full disk produced truncated output and exit status 0.
// stdio latches the error flag on the stream, so one std::ferror at the end
// catches every failed read or write on it, whenever it happened; and
// std::fclose reports separately, because the final buffer flush happens
// there and can fail on its own.
//
// The message cannot name the file: a unique_ptr deleter is stateless, and
// giving it a name would make FileHandle carry one for every stream. The
// point is that the run fails loudly rather than silently, which it did not
// before.
//
// This is also why the same deleter serves input and output. A read error
// was previously indistinguishable from end of input (Line_buffer sets
// at_end_ either way), so this reports that too.
//
// Calling fatal(), i.e. std::exit, from a destructor is safe here because
// every FileHandle is reached through a local in main() or in parse_fasta()
// and so is destroyed at normal scope exit. It would be undefined during
// std::exit's own static destruction, which is not a path any of them take.
// Defined in input_output.cpp, not here: it calls fatal(), and fatal.hpp
// has no include guard, so pulling it into a header that other headers
// reach would break any translation unit that arrives at it twice. Being
// out of line costs nothing -- it runs once per file, at close.
struct CloseFileHandle {
  auto operator()(std::FILE * file_handle) const -> void;
};

using FileHandle = std::unique_ptr<std::FILE, CloseFileHandle>;


// Both take the filename by reference rather than as a char const *: every
// caller holds a std::string, and a reference cannot be null, so the
// null-pointer contract these used to assert is now carried by the type.
auto fopen_input(std::string const & filename) -> FileHandle;
auto fopen_output(std::string const & filename) -> FileHandle;

#endif  // SWARM_UTILS_INPUT_OUTPUT_H
