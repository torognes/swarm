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

#include "swarm.h"
#include "utils/fatal.h"
#include <cstdint>  // uint64_t
#include <cstdio>  // FILE
#include <cstdlib>  // free, posix_memalign, realloc
#include <cstring>  // strcmp
#include <stdio.h>  // fdopen, ssize_t, getline
#include <stdlib.h> // posix_memalign
#include <unistd.h>  // dup, STDIN_FILENO, STDOUT_FILENO


static const char * progress_prompt;
static uint64_t progress_next;
static uint64_t progress_size;
static uint64_t progress_chunk;
constexpr std::size_t memalignment = 16;


auto progress_init(const char * prompt, const uint64_t size) -> void
{
  static constexpr uint64_t progress_granularity {200};
  progress_prompt = prompt;
  progress_size = size;
  progress_chunk = size < progress_granularity ?
    1 : size / progress_granularity;
  progress_next = 1;
  if (not opt_log.empty()) {
    std::fprintf(logfile, "%s", prompt);
  }
  else {
    std::fprintf(logfile, "%s %.0f%%", prompt, 0.0);
  }
}

auto progress_update(const uint64_t progress) -> void
{
  if (opt_log.empty() and (progress >= progress_next))
    {
      std::fprintf(logfile, "  \r%s %.0f%%", progress_prompt,
              100.0 * static_cast<double>(progress)
              / static_cast<double>(progress_size));
      progress_next = progress + progress_chunk;
      fflush(logfile);
    }
}

auto progress_done() -> void
{
  if (not opt_log.empty()) {
    std::fprintf(logfile, " %.0f%%\n", 100.0);
  }
  else {
    std::fprintf(logfile, "  \r%s %.0f%%\n", progress_prompt, 100.0);
  }
  fflush(logfile);
}


auto xmalloc(std::size_t size) -> void *
{
  if (size == 0) {
    size = 1;
  }
  void * t {nullptr};
#ifdef _WIN32
  t = _aligned_malloc(size, memalignment);
#else
  if (posix_memalign(& t, memalignment, size) != 0) {
    t = nullptr;
  }
#endif
  if (t == nullptr) {
    fatal(error_prefix, "Unable to allocate enough memory.");
  }
  return t;
}


auto xrealloc(void *ptr, std::size_t size) -> void *
{
  if (size == 0) {
    size = 1;
  }
#ifdef _WIN32
  void * t = _aligned_realloc(ptr, size, memalignment);
#else
  void * t = std::realloc(ptr, size);
#endif
  if (t == nullptr) {
    fatal(error_prefix, "Unable to allocate enough memory.");
  }
  return t;
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

auto fopen_input(const char * filename) -> std::FILE *
{
  /* open the input stream given by filename, but use stdin if name is - */
  std::FILE * input_stream = nullptr;

  if (std::strcmp(filename, "-") == 0) {
    const int file_descriptor = dup(STDIN_FILENO);
    input_stream = file_descriptor > 0 ? fdopen(file_descriptor, "rb") : nullptr;
  }
  else {
    input_stream = fopen(filename, "rb");  // refactoring: prefer std::fstream
  }

  return input_stream;
}

auto fopen_output(const char * filename) -> std::FILE *
{
  /* open the output stream given by filename, but use stdout if name is - */
  std::FILE * output_stream {nullptr};

  if (std::strcmp(filename, "-") == 0) {
    const int file_descriptor = dup(STDOUT_FILENO);
    output_stream = file_descriptor > 0 ? fdopen(file_descriptor, "w") : nullptr;
  }
  else {
    output_stream = fopen(filename, "w");  // refactoring: prefer std::fstream
  }

  return output_stream;
}

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
