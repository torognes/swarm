/*
    SWARM

    Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe

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
static const uint64_t progress_granularity = 200;
const size_t memalignment = 16;

void progress_init(const char * prompt, uint64_t size)
{
  progress_prompt = prompt;
  progress_size = size;
  progress_chunk = size < progress_granularity ?
    1 : size / progress_granularity;
  progress_next = 1;
  if (opt_log)
    fprintf(logfile, "%s", prompt);
  else
    fprintf(logfile, "%s %.0f%%", prompt, 0.0);
}

void progress_update(uint64_t progress)
{
  if ((!opt_log) && (progress >= progress_next))
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
  if (opt_log)
    fprintf(logfile, " %.0f%%\n", 100.0);
  else
    fprintf(logfile, "  \r%s %.0f%%\n", progress_prompt, 100.0);
  fflush(logfile);
}

int64_t gcd(int64_t a, int64_t b)
{
  if (b == 0)
  {
    return a;
  }
  else
  {
    return gcd(b, a % b);
  }
}

[[ noreturn ]] void fatal(const char * msg)
{
  fprintf(stderr, "\nError: %s\n", msg);
  exit(1);
}

void * xmalloc(size_t size)
{
  if (size == 0)
    size = 1;
  void * t = nullptr;
#ifdef _WIN32
  t = _aligned_malloc(size, memalignment);
#else
  if (posix_memalign(& t, memalignment, size))
    t = nullptr;
#endif
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  if (size == 0)
    size = 1;
#ifdef _WIN32
  void * t = _aligned_realloc(ptr, size, memalignment);
#else
  void * t = realloc(ptr, size);
#endif
  if (!t)
    fatal("Unable to reallocate enough memory.");
  return t;
}

void xfree(void * ptr)
{
  if (ptr)
    {
#ifdef _WIN32
      _aligned_free(ptr);
#else
      free(ptr);
#endif
    }
  else
    fatal("Trying to free a null pointer");
}

#if 0

/* never used */

uint64_t hash_fnv_1a_64(unsigned char * s, uint64_t n)
{
  const uint64_t fnv_offset = 14695981039346656037ULL;
  const uint64_t fnv_prime = 1099511628211; /* 2^40 - 435 */

  uint64_t hash = fnv_offset;

  for(uint64_t i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = (hash ^ c) * fnv_prime;
    }

  return hash;
}

unsigned int hash_fnv_1a_32(unsigned char * s, uint64_t n)
{
  const unsigned int fnv_offset = 2166136261;
  const unsigned int fnv_prime = 16777619;

  unsigned int hash = fnv_offset;

  for(uint64_t i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = (hash ^ c) * fnv_prime;
    }

  return hash;
}


uint64_t hash_djb2(unsigned char * s, uint64_t n)
{
  const uint64_t djb2_offset = 5381;

  uint64_t hash = djb2_offset;

  for(uint64_t i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = ((hash << 5) + hash) + c; /* hash = hash * 33 + c */
    }

  return hash;
}

uint64_t hash_djb2a(unsigned char * s, uint64_t n)
{
  const uint64_t djb2_offset = 5381;

  uint64_t hash = djb2_offset;

  for(uint64_t i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = ((hash << 5) + hash) ^ c; /* hash = hash * 33 ^ c */
    }

  return hash;
}

uint64_t hash64shift(uint64_t key)
{
  key = (~key) + (key << 21); // key = (key << 21) - key - 1;
  key = key ^ (key >> 24);
  key = (key + (key << 3)) + (key << 8); // key * 265
  key = key ^ (key >> 14);
  key = (key + (key << 2)) + (key << 4); // key * 21
  key = key ^ (key >> 28);
  key = key + (key << 31);
  return key;
}

uint64_t hash_xor64len(unsigned char * s, uint64_t n)
{
  uint64_t hash;

  hash = 8 * n;
  uint64_t * p = (uint64_t*) s;
  for(uint64_t i = 0; i < n/8; i++)
    hash ^= *p++;

  // Only the lowest (right-most) bits are used for indexing the hash table.
  // Make sure these bits are representative of all bits in the hash.

  hash = hash64shift(hash);

  return hash;
}

#endif

uint64_t hash_cityhash64(unsigned char * s, uint64_t n)
{
  return CityHash64(reinterpret_cast<const char *>(s), n);
}


FILE * fopen_input(const char * filename)
{
  /* open the input stream given by filename, but use stdin if name is - */
  if (strcmp(filename, "-") == 0)
    {
      int fd = dup(STDIN_FILENO);
      if (fd < 0)
        return nullptr;
      else
        return fdopen(fd, "rb");
    }
  else
    return fopen(filename, "rb");
}

FILE * fopen_output(const char * filename)
{
  /* open the output stream given by filename, but use stdout if name is - */
  if (strcmp(filename, "-") == 0)
    {
      int fd = dup(STDOUT_FILENO);
      if (fd < 0)
        return nullptr;
      else
        return fdopen(fd, "w");
    }
  else
    return fopen(filename, "w");
}
