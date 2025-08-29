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

#include "utils/fatal.h"
#include <cstdio>  // FILE // stdio.h: fdopen, ssize_t, getline
#include <cstdlib>  // free, posix_memalign, realloc
#include <cstring>  // strcmp


constexpr std::size_t memalignment = 16;


auto xmalloc(std::size_t size) -> void *
{
  if (size == 0) {
    size = 1;
  }
  void * memptr {nullptr};  // address of the allocated memory
#ifdef _WIN32
  memptr = _aligned_malloc(size, memalignment);
#else
  if (posix_memalign(& memptr, memalignment, size) != 0) {
    memptr = nullptr;
  }
#endif
  if (memptr == nullptr) {
    fatal(error_prefix, "Unable to allocate enough memory.");
  }
  return memptr;
}


auto xfree(void * ptr) -> void
{
  if (ptr != nullptr)
    {
#ifdef _WIN32
      _aligned_free(ptr);
#else
      std::free(ptr);
#endif
    }
  else {
    fatal(error_prefix, "Trying to free a null pointer.");
  }
}


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

  const std::size_t minsize = 2;
  const std::size_t maxsize = SIZE_MAX / 2;

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
      *linep = (char *) malloc(*linecapp);  // refactor C++11: *linep = new char[*linecapp];
      if (*linep == nullptr)
        return -1;
    }

  char * p = *linep;            // pointer to where to put next char
  char * e = p + *linecapp - 1; // pointer to last byte in buffer
  std::size_t len = 0;
  *p = 0;

  while (1)
    {
      while (p < e)
        {
          int c = getc(stream);
          switch(c)
            {
            case -1:
              if (feof(stream))
                {
                  // EOF, add NUL
                  *p = 0;
                  len = p - *linep;
                  if (len > 0)
                    return len;
                  else
                    return -1;
                }
              else
                {
                  // Error
                  return -1;
                }

            case '\n':
              // Newline
              *p = c;
              ++p;
              *p = 0;
              return p - *linep;

            default:
              // Ordinary character, including NUL
              *p = c;
              ++p;
              break;
            }
        }

      // Increase buffer size

      if (*linecapp >= maxsize)
        {
          errno = EOVERFLOW;
          return -1;
        }

      std::size_t newlinecap = minsize;
      while ((newlinecap <= *linecapp) and (newlinecap < maxsize))
        newlinecap *= 2;

      if (newlinecap > maxsize)
        {
          errno = EOVERFLOW;
          return -1;
        }

      char * newlinep = (char *) std::realloc(*linep, newlinecap);
      if (newlinep == nullptr)
        {
          // Memory allocation error
          return -1;
        }

      len = p - *linep;
      *linep = newlinep;
      *linecapp = newlinecap;
      p = newlinep + len;
      e = p + *linecapp - 1;
    }
#endif
}
