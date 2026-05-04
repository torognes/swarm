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

#include "fatal.h"
#include "xgetline.h"
#include <algorithm>  // std::min
#include <cstdio>  // FILE // stdio.h: fdopen, ssize_t, getline
#include <cstdlib>  // malloc, realloc, free
#include <iterator>  // std::next
#include <string>  // std::char_traits


// refactoring: std::getline(input, str) -> input
auto xgetline(char ** linep, std::size_t * linecapp, std::FILE * stream) -> ssize_t
{
#ifndef _WIN32

  return getline(linep, linecapp, stream);

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

  static constexpr std::size_t minsize {2};
  static constexpr std::size_t maxsize {SIZE_MAX / 2};
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
      *linecapp = minsize;
      *linep = static_cast<char *>(std::malloc(*linecapp));
      if (*linep == nullptr) {
        return -1;
      }
    }

  auto * p = *linep;                  // pointer to where to put next char
  auto const * e = std::next(p, static_cast<std::ptrdiff_t>(*linecapp - 1));  // pointer to last byte in buffer
  *p = '\0';

  while (true)
    {
      while (p < e)
        {
          auto const c = std::getc(stream);
          if (c == eof_value)
            {
              if (std::feof(stream) != 0)
                {
                  // EOF, add NUL
                  *p = '\0';
                  auto const len = static_cast<std::size_t>(p - *linep);
                  if (len > 0) {
                    return static_cast<ssize_t>(len);
                  }
                  return -1;
                }
              // Error
              return -1;
            }
          if (c == '\n')
            {
              // Newline
              *p = static_cast<char>(c);
              ++p;
              *p = '\0';
              return p - *linep;
            }
          // Ordinary character, including NUL
          *p = static_cast<char>(c);
          ++p;
        }

      // Increase buffer size

      if (*linecapp >= maxsize)
        {
          errno = EOVERFLOW;
          return -1;
        }

      auto const newlinecap = std::min(*linecapp * 2, maxsize);

      auto * const newlinep = static_cast<char *>(std::realloc(*linep, newlinecap));
      if (newlinep == nullptr)
        {
          // Memory allocation error
          return -1;
        }

      auto const len = static_cast<std::size_t>(p - *linep);
      *linep = newlinep;
      *linecapp = newlinecap;
      p = std::next(newlinep, static_cast<std::ptrdiff_t>(len));
      e = std::next(p, static_cast<std::ptrdiff_t>(*linecapp - 1));
    }
#endif
}


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
  auto const linelen = xgetline(&data_, &capacity_, stream);
  if (linelen < 0) {
    *data_ = '\0';
    return;
  }
  filepos += static_cast<unsigned long int>(linelen);
}
