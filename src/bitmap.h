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

class Bitmap
{
 private:
  size_t size;      /* size in bits */
  unsigned char * data;    /* the actual bitmap */
  
 public:

  explicit Bitmap(size_t _size)
  {
    size = _size;
    data = (unsigned char *) xmalloc((size+7)/8);
  }
  
  ~Bitmap()
  {
    if (data)
      free(data);
  }
  
  bool get(size_t x)
  {
    return (data[x >> 3] >> (x & 7)) & 1;
  }

  void reset_all()
  {
    memset(data, 0, (size+7)/8);
  }
  
  void set_all()
  {
    memset(data, 255, (size+7)/8);
  }
  
  void reset(size_t x)
  {
    //    data[x >> 3] &= ~ (1 << (x & 7));
    __sync_fetch_and_and(data + (x >> 3), ~(1 << (x & 7)));
  }
  
  void set(size_t x)
  {
    //    data[x >> 3] |= 1 << (x & 7);
    __sync_fetch_and_or(data + (x >> 3), 1 << (x & 7));
  }
  
  void flip(size_t x)
  {
    //    data[x >> 3] ^= 1 << (x & 7);
    __sync_fetch_and_xor(data + (x >> 3), 1 << (x & 7));
  }
};
