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

#ifndef SWARM_UTILS_ELEMENT_ORDER_H
#define SWARM_UTILS_ELEMENT_ORDER_H


// How View and Span order their elements.
//
// The primary template defers to the element type's own operator<, which is
// what a bare std::lexicographical_compare would have done anyway.
//
// The specialization for char is the reason this trait exists: it orders bytes
// as unsigned char. That is what std::strcmp does, and what std::string does
// through std::char_traits<char>::lt -- but it is *not* what comparing char
// with operator< does, because char is signed on x86-64 and on the Windows
// target while ARM and PowerPC Linux default it to unsigned. Without the
// specialization, View<char>{} < View<char>{} would order any byte with its
// high bit set differently from std::string, differently from std::strcmp, and
// differently from one architecture to the next.
//
// In swarm that divergence was not hypothetical, it was observed. Sequence
// labels are arbitrary bytes and routinely carry UTF-8 or Latin-1 (an accented
// author name, a locality, a non-ASCII sample tag), and the label is the
// tie-break of three sort comparators: compare_entries (db.cpp), which fixes
// the cluster order written to -o, -i, -s and -u, and compare_seeds
// (algo_output.cpp, d > 1) and compare_mass_and_headers
// (algod1_output.cpp, d = 1), which fix the order of -w. So
//
//   >Frad_1 >Frzd_1 >Fréd_1
//
// came out as Fréd, Frad, Frzd on x86-64 and Windows, and as Frad, Frzd, Fréd
// on ARM and PowerPC: the same input, the same swarm version, two different
// outputs. The second of those is the strcmp order, which is also what swarm
// printed before the label comparison moved from std::strcmp to View::operator<.
//
// Only char is specialized. signed char and unsigned char say what they mean
// and are left to the primary template, the same line std::char_traits draws.

template <typename Type>
struct element_order {
  static constexpr auto less(Type const & lhs, Type const & rhs) -> bool {
    return lhs < rhs;
  }

  // Three-way, with the sign convention of std::strcmp: negative if lhs sorts
  // first, positive if rhs does, zero if the two are equivalent. Expressed with
  // two less() calls because the primary template can assume nothing beyond a
  // strict weak ordering.
  // C++14 refactoring: constexpr (C++11 allows a single return statement
  // only, which neither compare() can express without losing its shape)
  static auto compare(Type const & lhs, Type const & rhs) -> int {
    if (less(lhs, rhs)) { return -1; }
    // the swapped operands are the point: this asks the same question the
    // other way round, which is how a three-way result comes out of a
    // two-way predicate
    // NOLINTNEXTLINE(readability-suspicious-call-argument)
    if (less(rhs, lhs)) { return +1; }
    return 0;
  }
};

template <>
struct element_order<char> {
  static constexpr auto less(char const lhs, char const rhs) -> bool {
    return static_cast<unsigned char>(lhs) < static_cast<unsigned char>(rhs);
  }

  // One byte comparison rather than the primary template's two: char is
  // totally ordered once it is read as unsigned char.
  static auto compare(char const lhs, char const rhs) -> int {
    auto const lhs_byte = static_cast<unsigned char>(lhs);
    auto const rhs_byte = static_cast<unsigned char>(rhs);
    if (lhs_byte == rhs_byte) { return 0; }
    return (lhs_byte < rhs_byte) ? -1 : +1;
  }
};


// Functor wrapper, for the std algorithms that take a comparison object.
template <typename Type>
struct element_less {
  auto operator()(Type const & lhs, Type const & rhs) const -> bool {
    return element_order<Type>::less(lhs, rhs);
  }
};

#endif // SWARM_UTILS_ELEMENT_ORDER_H
