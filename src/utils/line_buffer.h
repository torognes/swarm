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

#include <cstdint>  // uint64_t
#include <cstdio>  // size_t


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
  // by the number of bytes consumed. On read failure, the buffer is
  // left empty (first byte set to '\0'); callers can use empty() as
  // the end-of-input sentinel.
  auto read_next(std::FILE * stream, uint64_t & filepos) -> void;

  auto data()       const noexcept -> char const * { return data_; }
  auto empty()      const noexcept -> bool         { return *data_ == '\0'; }
  auto peek_first() const noexcept -> char         { return *data_; }

private:
  // read_next() is the only writer to data_; it goes through the
  // field directly. Keeping data() public-const-only means external
  // callers cannot get a writable pointer into the buffer.
  char *      data_     {nullptr};
  std::size_t capacity_ {0};
};
