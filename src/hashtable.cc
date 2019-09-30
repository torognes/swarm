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

#define HASHFILLPCT 70

uint64_t hash_mask;
unsigned char * hash_occupied = nullptr;
uint64_t * hash_values = nullptr;
unsigned int * hash_data = nullptr;
uint64_t hash_tablesize = 0;

void hash_zap()
{
  memset(hash_occupied, 0, (hash_tablesize + 63) / 8);
}

void hash_alloc(uint64_t amplicons)
{
  hash_tablesize = 1;
  while (100 * amplicons > HASHFILLPCT * hash_tablesize)
    hash_tablesize <<= 1;
  hash_mask = hash_tablesize - 1;

  hash_occupied =
    static_cast<unsigned char *>(xmalloc((hash_tablesize + 63) / 8));
  hash_zap();

  hash_values =
    static_cast<uint64_t *>(xmalloc(hash_tablesize * sizeof(uint64_t)));

  hash_data = static_cast<unsigned int *>
    (xmalloc(hash_tablesize * sizeof(unsigned int)));
}

void hash_free()
{
  xfree(hash_occupied);
  xfree(hash_values);
  xfree(hash_data);
}
