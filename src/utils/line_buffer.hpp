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

#ifndef SWARM_UTILS_LINE_BUFFER_H
#define SWARM_UTILS_LINE_BUFFER_H


#include "view.hpp"  // View<char>
#include <cstddef>  // std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // std::FILE


// RAII wrapper for the line buffer used by read_next(). POSIX
// getline() owns the buffer's lifetime: it may std::realloc() it on
// long lines, so the storage must come from std::malloc and the
// destructor must call std::free. std::vector<char> or new[]/delete[]
// would create an allocator mismatch and undefined behavior.
class Line_buffer {
public:
  explicit Line_buffer(std::size_t initial);

  ~Line_buffer() noexcept;

  auto release() noexcept -> void;

  Line_buffer(Line_buffer const &)                     = delete;
  auto operator=(Line_buffer const &) -> Line_buffer & = delete;
  Line_buffer(Line_buffer &&)                          = delete;
  auto operator=(Line_buffer &&)      -> Line_buffer & = delete;

  // Read one line from `stream` into the buffer and bump `filepos`
  // by the number of bytes consumed. On read failure the buffer is left
  // empty and at_end() becomes true.
  auto read_next(std::FILE * stream, uint64_t & filepos) -> void;

  auto data()       const noexcept -> char const * { return data_; }
  auto peek_first() const noexcept -> char         { return *data_; }

  // True once the stream is exhausted. Deliberately distinct from "this
  // line has no content": a line that begins with a null byte has zero
  // visible length (see read_next) without the input being over, and
  // conflating the two made swarm discard the rest of the file.
  auto at_end() const noexcept -> bool { return at_end_; }

  // The line as read, without the terminating '\0'. read_next() stores
  // the length instead of discarding it, so consumers no longer have to
  // re-derive it by scanning for the sentinel.
  auto view() const noexcept -> View<char> { return View<char>{data_, length_}; }

private:
  // read_next() is the only writer to data_; it goes through the
  // field directly. Keeping data() public-const-only means external
  // callers cannot get a writable pointer into the buffer.
  char *      data_     {nullptr};
  std::size_t capacity_ {0};
  std::size_t length_   {0};
  bool        at_end_   {false};
};

#endif  // SWARM_UTILS_LINE_BUFFER_H
