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
#include "util.h"
#include "utils/hashtable_size.h"


uint64_t hash_mask {0};
unsigned char * hash_occupied {nullptr};
uint64_t * hash_values {nullptr};
unsigned int * hash_data {nullptr};
uint64_t hashtablesize {0};

void hash_zap()
{
  constexpr int padding {63};  // make sure our final value is >= 64 / 8
  constexpr int convert_to_bytes {8};
  memset(hash_occupied, 0, (hashtablesize + padding) / convert_to_bytes);
}

void hash_alloc(const uint64_t amplicons)
{
  constexpr int padding {63};  // make sure our final value is >= 64 / 8
  constexpr int convert_to_bytes {8};

  hashtablesize = compute_hashtable_size(amplicons);
  hash_mask = hashtablesize - 1;

  hash_occupied =
    static_cast<unsigned char *>(xmalloc((hashtablesize + padding) / convert_to_bytes));

  hash_zap();

  hash_values =
    static_cast<uint64_t *>(xmalloc(hashtablesize * sizeof(uint64_t)));

  hash_data = static_cast<unsigned int *>
    (xmalloc(hashtablesize * sizeof(unsigned int)));
}

void hash_free()
{
  xfree(hash_occupied);
  xfree(hash_values);
  xfree(hash_data);
}
