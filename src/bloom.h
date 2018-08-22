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

class BloomFilter
{
private:

  Bitmap bitmap;
  size_t m; /* total number of bits in bitmap */
  int k; /* number of hash functions */

public:

  BloomFilter(unsigned long _m, int _k) : bitmap(_m)
  {
    bitmap.reset_all();
    m = _m;
    k = _k;
  }

  bool get(const char * buf, size_t len)
  {
    uint128 hash = CityHash128(buf, len);
    uint64 h0 = Uint128Low64(hash);
    uint64 h1 = Uint128High64(hash);

    for(int i=0; i<k; i++)
      {
        uint64 h = h0 ^ (i * h1);
        if (! bitmap.get(h % m))
          return false;
      }
    return true;
  }

  void set(const char * buf, size_t len)
  {
    uint128 hash = CityHash128(buf, len);
    uint64 h0 = Uint128Low64(hash);
    uint64 h1 = Uint128High64(hash);

    for(int i=0; i<k; i++)
      {
        uint64 h = h0 ^ (i * h1);
        bitmap.set(h % m);
      }
  }
};
