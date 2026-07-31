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

#ifndef SWARM_UTILS_CHAIN_RANGE_H
#define SWARM_UTILS_CHAIN_RANGE_H

#include <cstddef>  // std::ptrdiff_t
#include <iterator>  // std::input_iterator_tag


// A range over an intrusive singly-linked chain held in an index table:
// each id names the next id, until a sentinel ends the chain. swarm keeps
// two of these -- the amplicons of one cluster (ampinfo_s::next, ending at
// no_swarm) and the identical copies of a dereplicated sequence
// (nextseqtab[id], ending at 0) -- and walked both by hand at eight call
// sites, each repeating the "id = table[id].next" step and the sentinel.
//
// Next is a policy type supplying
//   static auto next(Table const & table, Id id) -> Id
// so that a chain threaded through a field (ampinfo_s::next) and one
// threaded through a flat table (nextseqtab) share this iteration without
// this header knowing either type.
//
// The iterator yields ids by value, which in C++11 makes it an input
// iterator rather than a forward one (its reference type is not a real
// reference, the same reason std::vector<bool>'s iterator is not a forward
// iterator). That is enough for the range-for loops these chains are used
// in; it is deliberately not enough for std::search and friends.
template <typename Table, typename Next, typename Id = unsigned int>
class Chain_range {
public:
  Chain_range(Table const & table, Id const first, Id const sentinel) noexcept
    : table_ {&table},
      first_ {first},
      sentinel_ {sentinel} {
  }

  class const_iterator {
  public:
    using iterator_category = std::input_iterator_tag;
    using value_type        = Id;
    using difference_type   = std::ptrdiff_t;
    using pointer           = Id const *;
    using reference         = Id;

    const_iterator(Table const * const table, Id const current) noexcept
      : table_ {table},
        current_ {current} {
    }

    auto operator*() const noexcept -> Id { return current_; }

    auto operator++() noexcept -> const_iterator & {
      current_ = Next::next(*table_, current_);
      return *this;
    }

    // Compares the position in the chain only: two iterators into the same
    // chain agree on the table, and end() is the sentinel position.
    auto operator==(const_iterator const & other) const noexcept -> bool {
      return current_ == other.current_;
    }
    auto operator!=(const_iterator const & other) const noexcept -> bool {
      return not (*this == other);
    }

  private:
    Table const * table_;
    Id            current_;
  };

  auto begin() const noexcept -> const_iterator {
    return const_iterator{table_, first_};
  }
  auto end() const noexcept -> const_iterator {
    return const_iterator{table_, sentinel_};
  }

private:
  Table const * table_;
  Id            first_;
  Id            sentinel_;
};

#endif  // SWARM_UTILS_CHAIN_RANGE_H
