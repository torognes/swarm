/*
    SWARM

    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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

static const char * progress_prompt;
static uint64_t progress_next;
static uint64_t progress_size;
static uint64_t progress_chunk;
const size_t memalignment = 16;

void progress_init(const char * prompt, uint64_t size)
{
  constexpr uint64_t progress_granularity {200};
  progress_prompt = prompt;
  progress_size = size;
  progress_chunk = size < progress_granularity ?
    1 : size / progress_granularity;
  progress_next = 1;
  if (! opt_log.empty()) {
    fprintf(logfile, "%s", prompt);
  }
  else {
    fprintf(logfile, "%s %.0f%%", prompt, 0.0);
  }
}

void progress_update(uint64_t progress)
{
  if (opt_log.empty() && (progress >= progress_next))
    {
      fprintf(logfile, "  \r%s %.0f%%", progress_prompt,
              100.0 * static_cast<double>(progress)
              / static_cast<double>(progress_size));
      progress_next = progress + progress_chunk;
      fflush(logfile);
    }
}

void progress_done()
{
  if (! opt_log.empty()) {
    fprintf(logfile, " %.0f%%\n", 100.0);
  }
  else {
    fprintf(logfile, "  \r%s %.0f%%\n", progress_prompt, 100.0);
  }
  fflush(logfile);
}

// std::gcd() in C++17
auto gcd(int64_t a, int64_t b) -> int64_t
{
  return b == 0 ? a : gcd(b, a % b);
}

[[ noreturn ]] void fatal(const char * msg)
{
  fprintf(stderr, "\nError: %s\n", msg);
  exit(1);
}

auto xmalloc(size_t size) -> void *
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
    fatal("Unable to allocate enough memory.");
  }
  return t;
}

auto xrealloc(void *ptr, size_t size) -> void *
{
  if (size == 0) {
    size = 1;
  }
#ifdef _WIN32
  void * t = _aligned_realloc(ptr, size, memalignment);
#else
  void * t = realloc(ptr, size);
#endif
  if (t == nullptr) {
    fatal("Unable to reallocate enough memory.");
  }
  return t;
}

void xfree(void * ptr)
{
  if (ptr != nullptr)
    {
#ifdef _WIN32
      _aligned_free(ptr);
#else
      free(ptr);
#endif
    }
  else {
    fatal("Trying to free a null pointer");
  }
}

auto fopen_input(const char * filename) -> std::FILE *
{
  /* open the input stream given by filename, but use stdin if name is - */
  std::FILE * input_stream = nullptr;

  if (strcmp(filename, "-") == 0) {
    int fd = dup(STDIN_FILENO);
    input_stream = fd > 0 ? fdopen(fd, "rb") : nullptr;
  }
  else {
    input_stream = fopen(filename, "rb");
  }

  return input_stream;
}

auto fopen_output(const char * filename) -> std::FILE *
{
  /* open the output stream given by filename, but use stdout if name is - */
  std::FILE * output_stream {nullptr};

  if (strcmp(filename, "-") == 0) {
    int fd = dup(STDOUT_FILENO);
    output_stream = fd > 0 ? fdopen(fd, "w") : nullptr;
  }
  else {
    output_stream = fopen(filename, "w");
  }

  return output_stream;
}

ssize_t xgetline(char ** linep, size_t * linecapp, FILE * stream)
{
  /*
     Replacement for the getline function.
     May be used on Windows and other non-POSIX systems.
     Dynamic buffer expansion while reading input using fgets.
  */

  const size_t minsize = 2;
  const size_t maxsize = (size_t)(SSIZE_MAX) + 1;

  /* Error if linep or linecapp pointers are null */
  if ((linep == nullptr) || (linecapp == nullptr))
    {
      return -1;
    }

  if (*linep == nullptr)
    {
      /* allocate a default buffer if linep is a null pointer */
      *linecapp = minsize;
      *linep = (char *) xmalloc(*linecapp);
    }
  else if (*linecapp < minsize)
    {
      /* reallocate a default buffer if provided buffer is too small */
      *linecapp = minsize;
      *linep = (char *) xrealloc(*linep, *linecapp);
    }

  /* repeatedly call fgets until entire line has been read */

  char * p = *linep;      // pointer to where to read next
  size_t len = 0;         // total number of chars read
  size_t rem = *linecapp; // remaining capacity in buffer

  while (1)
    {
      /* size = amount to read using fgets */
      int size;
      if (rem <= INT_MAX)
        size = rem;
      else
        size = INT_MAX;

      /* read using fgets */
      char * s = fgets(p, size, stream);

      /* EOF before any characters read? or Error? */
      if ((s == nullptr) || ferror(stream))
        {
          return -1;
        }

      /* no error, at least 1 char read */

      /* determine length of string just read */
      uint64_t slen = strlen(s);
      p += slen;
      len += slen;
      rem -= slen;

      /* check if last character was newline */
      /* or if we have reached eof */
      if ((s[slen-1] == '\n') || feof(stream))
        {
          /* Success: we are done, return length */
          return len;
        }
      else
        {
          /* Need to read more */

          if (len > maxsize)
            {
              /* Error: Delimiter not found in first SSIZE_MAX chars. */
              return -1;
            }

          /* Expand buffer? Needs space for at least 1 char + NUL */
          if (rem < minsize)
            {
              /* Already too big */
              if (*linecapp >= maxsize)
                return -1;

              /* Increase buffer size exponentially */
              size_t newlinecap = minsize;
              while ((newlinecap <= *linecapp) && (newlinecap < maxsize))
                newlinecap *= 2;

              if (newlinecap > maxsize)
                return -1;

              /* Expand buffer using realloc */
              char * newlinep = (char *) realloc(*linep, newlinecap);
              if (newlinep == nullptr)
                {
                  /* memory allocation error */
                  return -1;
                }

              /* updates */
              size_t expand = newlinecap - *linecapp;
              *linep = newlinep;
              *linecapp = newlinecap;
              p = newlinep + len;
              rem += expand;
            }

          /* read more */
        }
    }
}
