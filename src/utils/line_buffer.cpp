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

#include "fatal.hpp"
#include "line_buffer.hpp"
#include <algorithm>  // std::find // _WIN32: std::min
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // uint64_t
#include <cstdio>  // FILE // stdio.h: fdopen, ssize_t, getline
#include <cstdlib>  // malloc, realloc, free
#include <iterator>  // std::distance, std::next

#ifdef _WIN32
#include <cerrno>  // errno, EINVAL, EOVERFLOW
#include <limits>  // std::numeric_limits
#include <string>  // std::char_traits
#endif


namespace {

// refactoring: std::getline(input, str) -> input
//
// The non-POSIX fallback below is required because POSIX getline()
// is not part of the C++ standard library and is not provided by
// MinGW/Windows toolchains. Swarm supports Windows builds (see the
// mingw target in the Makefile), so a portable replacement must be
// available when _WIN32 is defined.
auto read_one_line(char ** linep, std::size_t * linecapp, std::FILE * stream) -> ssize_t
{
#ifndef _WIN32

  return ::getline(linep, linecapp, stream);

#else

  /*
     Replacement for the POSIX getline function.
     May be used on Windows and other non-POSIX systems.
     Dynamic buffer expansion while reading input.
     Considerably slower since it calls getc repeatedly.

     Using fgets is much faster but does not work properly as
     it cannot handle NUL characters in the string correctly,
     which is important for correct counting of characters and file size.
  */

  static constexpr std::size_t min_capacity {2};
  static constexpr std::size_t max_capacity {std::numeric_limits<std::size_t>::max() / 2};
  static constexpr auto eof_value = std::char_traits<char>::eof();

  /* Error if linep or linecapp pointers are null */
  if ((linep == nullptr) or (linecapp == nullptr))
    {
      errno = EINVAL;
      return -1;
    }

  if (*linep == nullptr)
    {
      /* allocate a default buffer if linep is a null pointer */
      *linecapp = min_capacity;
      *linep = static_cast<char *>(std::malloc(*linecapp));
      if (*linep == nullptr) {
        return -1;
      }
    }

  auto * cursor = *linep;                  // pointer to where to put next char
  auto const * end_of_buffer = std::next(cursor, static_cast<std::ptrdiff_t>(*linecapp - 1));  // pointer to last byte in buffer
  *cursor = '\0';

  while (true)
    {
      while (cursor < end_of_buffer)
        {
          auto const character = std::getc(stream);
          if (character == eof_value)
            {
              if (std::feof(stream) != 0)
                {
                  // EOF, add NUL
                  *cursor = '\0';
                  auto const length = static_cast<std::size_t>(cursor - *linep);
                  if (length > 0) {
                    return static_cast<ssize_t>(length);
                  }
                  return -1;
                }
              // Error
              return -1;
            }
          if (character == '\n')
            {
              // Newline
              *cursor = static_cast<char>(character);
              ++cursor;
              *cursor = '\0';
              return cursor - *linep;
            }
          // Ordinary character, including NUL
          *cursor = static_cast<char>(character);
          ++cursor;
        }

      // Increase buffer size

      if (*linecapp >= max_capacity)
        {
          errno = EOVERFLOW;
          return -1;
        }

      auto const new_capacity = std::min(*linecapp * 2, max_capacity);

      auto * const new_buffer = static_cast<char *>(std::realloc(*linep, new_capacity));
      if (new_buffer == nullptr)
        {
          // Memory allocation error
          return -1;
        }

      auto const length = static_cast<std::size_t>(cursor - *linep);
      *linep = new_buffer;
      *linecapp = new_capacity;
      cursor = std::next(new_buffer, static_cast<std::ptrdiff_t>(length));
      end_of_buffer = std::next(new_buffer, static_cast<std::ptrdiff_t>(*linecapp - 1));
    }
#endif
}

}  // namespace


Line_buffer::Line_buffer(std::size_t const initial)
  : data_{static_cast<char *>(std::malloc(initial))}, capacity_{initial}
{
  if (data_ == nullptr) {
    fatal("Unable to allocate enough memory.");
  }
}


Line_buffer::~Line_buffer() noexcept { release(); }


auto Line_buffer::release() noexcept -> void {
  if (data_ != nullptr) {
    std::free(data_);
    data_ = nullptr;
    capacity_ = 0;
  }
}


auto Line_buffer::read_next(std::FILE * stream, uint64_t & filepos) -> void
{
  auto const linelen = read_one_line(&data_, &capacity_, stream);
  if (linelen < 0) {
    *data_ = '\0';
    length_ = 0;
    return;
  }
  filepos += static_cast<unsigned long int>(linelen);

  // read_one_line() stores a NUL byte read from the input verbatim (see
  // the _WIN32 branch above, and getline() likewise), so a line may
  // contain one. The published length stops at the first such byte,
  // which is where every consumer stopped when it walked to the '\0'
  // sentinel itself: this keeps the truncation swarm has always
  // performed rather than making the bytes after a NUL newly visible.
  // filepos still counts the whole line, since it tracks the file
  // position and not what was read out of the buffer.
  //
  // Whether such a line should be rejected outright instead of silently
  // truncated is a separate, open question: an existing test asserts
  // that it is accepted (test_input.sh, "ascii character 0 is allowed
  // in sequences"), and the manual page says nothing about NUL bytes.
  auto const * const line_begin = data_;
  auto const * const line_end = std::next(line_begin, static_cast<std::ptrdiff_t>(linelen));
  auto const * const first_nul = std::find(line_begin, line_end, '\0');
  length_ = static_cast<std::size_t>(std::distance(line_begin, first_nul));
}
