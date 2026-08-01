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

#include <array>
#include <cstdint> // uint64_t, std::uint8_t
#include <vector>


struct Sequence;  // defined in db.hpp


class Zobrist {
public:
  explicit Zobrist(unsigned int n);

  auto value(unsigned int pos, unsigned char offset) const -> uint64_t;

  auto hash(Sequence const & seq) const -> uint64_t;
  auto hash_delete_first(Sequence const & seq) const -> uint64_t;
  auto hash_insert_first(Sequence const & seq) const -> uint64_t;

private:
  enum struct First_base_op : std::uint8_t { remove, insert_gap };
  auto hash_first_shifted(Sequence const & seq, First_base_op operation) const -> uint64_t;

  auto fill_rng_table(unsigned int zobrist_len) -> void;
  auto fill_rng_byte_table(unsigned int zobrist_len) -> void;

  static constexpr auto nt_per_byte = 4U;     // 4 nucleotides packed per encoded byte
  static constexpr auto byte_range = 256U;    // 8-bit byte values: 256 possibilities

  // one row of tab_byte_base_v_ below: named so that hash() can spell the
  // element type of the range it zips the encoded bytes against
  using Byte_row = std::array<uint64_t, byte_range>;

  // tab_base_v_[pos]      : 4 RNG values, one per nucleotide A/C/G/T
  // tab_byte_base_v_[bpos]: 256 precomputed XOR-folds, one per possible byte value,
  //                         where bpos is the byte position in the encoded buffer
  std::vector<std::array<uint64_t, nt_per_byte>>  tab_base_v_;
  std::vector<Byte_row>                           tab_byte_base_v_;
};


