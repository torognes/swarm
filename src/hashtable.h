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

extern uint64_t hash_mask;
extern unsigned char * hash_occupied;
extern uint64_t * hash_values;
extern int * hash_data;
extern uint64_t hash_tablesize;

inline uint64_t hash_get_tablesize()
{
  return hash_tablesize;
}

inline uint64_t hash_getindex(uint64_t hash)
{
  // Shift bits right to get independence from the simple Bloom filter hash
  hash = hash >> 32;
  return hash & hash_mask;
}

inline unsigned int hash_getnextindex(unsigned int j)
{
  return (j+1) & hash_mask;
}

inline void hash_set_occupied(unsigned int j)
{
  hash_occupied[j >> 3] |= (1 << (j & 7));
}

inline int hash_is_occupied(unsigned int j)
{
  return hash_occupied[j >> 3] & (1 << (j & 7));
}

inline void hash_set_value(unsigned int j, uint64_t hash)
{
  hash_values[j] = hash;
}

inline int hash_compare_value(unsigned int j, uint64_t hash)
{
  return (hash_values[j] == hash);
}

inline int hash_get_data(unsigned int j)
{
  return hash_data[j];
}

inline void hash_set_data(unsigned int j, int x)
{
  hash_data[j] = x;
}

void hash_zap();

void hash_alloc(uint64_t amplicons);

void hash_free();
