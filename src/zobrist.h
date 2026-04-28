/*
    SWARM

    Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

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

#include <cstdint> // uint64_t
#include <vector>


class Zobrist {
public:
  explicit Zobrist(unsigned int n);

  auto value(unsigned int pos, unsigned char offset) const -> uint64_t;
  auto hash(char const * seq, unsigned int len) const -> uint64_t;
  auto hash_delete_first(char const * seq, unsigned int len) const -> uint64_t;
  auto hash_insert_first(char const * seq, unsigned int len) const -> uint64_t;

private:
  std::vector<uint64_t> tab_base_v_;
  std::vector<uint64_t> tab_byte_base_v_;
};


