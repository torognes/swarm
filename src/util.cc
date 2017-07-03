/*
    SWARM

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

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
static unsigned long progress_next;
static unsigned long progress_size;
static unsigned long progress_chunk;
static const unsigned long progress_granularity = 200;

void progress_init(const char * prompt, unsigned long size)
{
  progress_prompt = prompt;
  progress_size = size;
  progress_chunk = size < progress_granularity ? 
    1 : size / progress_granularity;
  progress_next = 0;
  if (opt_log)
    fprintf(logfile, "%s", prompt);
  else
    fprintf(logfile, "%s %.0f%%", prompt, 0.0);
}

void progress_update(unsigned long progress)
{
  if ((!opt_log) && (progress >= progress_next))
    {
      fprintf(logfile, "  \r%s %.0f%%", progress_prompt,
              100.0 * progress / progress_size);
      progress_next = progress + progress_chunk;
    }
}

void progress_done()
{
  if (opt_log)
    fprintf(logfile, " %.0f%%\n", 100.0);
  else
    fprintf(logfile, "  \r%s %.0f%%\n", progress_prompt, 100.0);
}

long gcd(long a, long b)
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

void fatal(const char * msg)
{
  fprintf(stderr, "\nError: %s\n", msg);
  exit(1);
}

void fatal(const char * format, const char * message)
{
  fprintf(stderr, "\n");
  fprintf(stderr, format, message);
  fprintf(stderr, "\n");
  exit(1);
}

void * xmalloc(size_t size)
{
  const size_t alignment = 16;
  void * t = NULL;
  if (posix_memalign(& t, alignment, size))
    fatal("Unable to allocate enough memory.");
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

void * xrealloc(void *ptr, size_t size)
{
  void * t = realloc(ptr, size);
  if (!t)
    fatal("Unable to allocate enough memory.");
  return t;
}

unsigned long hash_fnv_1a_64(unsigned char * s, unsigned long n)
{
  const unsigned long fnv_offset = 14695981039346656037UL;
  const unsigned long fnv_prime = 1099511628211; /* 2^40 - 435 */

  unsigned long hash = fnv_offset;

  for(unsigned long i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = (hash ^ c) * fnv_prime;
    }

  return hash;
}

unsigned int hash_fnv_1a_32(unsigned char * s, unsigned long n)
{
  const unsigned int fnv_offset = 2166136261;
  const unsigned int fnv_prime = 16777619;

  unsigned int hash = fnv_offset;

  for(unsigned long i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = (hash ^ c) * fnv_prime;
    }

  return hash;
}

unsigned long hash_djb2(unsigned char * s, unsigned long n)
{
  const unsigned long djb2_offset = 5381;

  unsigned long hash = djb2_offset;

  for(unsigned long i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = ((hash << 5) + hash) + c; /* hash = hash * 33 + c */
    }

  return hash;
}

unsigned long hash_djb2a(unsigned char * s, unsigned long n)
{
  const unsigned long djb2_offset = 5381;

  unsigned long hash = djb2_offset;

  for(unsigned long i = 0; i < n; i++)
    {
      unsigned char c = *s++;
      hash = ((hash << 5) + hash) ^ c; /* hash = hash * 33 ^ c */
    }

  return hash;
}

unsigned long hash_cityhash64(unsigned char * s, unsigned long n)
{
  return CityHash64((const char*)s, n);
}
