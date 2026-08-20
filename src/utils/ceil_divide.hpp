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

#ifndef SWARM_UTILS_CEIL_DIVIDE_H
#define SWARM_UTILS_CEIL_DIVIDE_H

#include <cassert>
#include <cstdint>  // uint64_t
#include <type_traits>

#ifndef NDEBUG
#include <limits>
#endif


// "Ceiling division": smallest integer >= numerator / denominator.
//
// idiom:
//   (numerator + (denominator - 1)) / denominator
//
// example: round 11 / 8 up to 2
//
//   11 + 7  = 18
//   18 / 8  = 2  (integer division truncates the .25 remainder)
//
// note: if numerator > max - (denominator - 1), the intermediate
//       sum overflows => Error!
//
// The sum-first idiom is kept for its cost, not only its brevity: where the
// numerator arrives as a sum, as it does in qgram_diff (two psadbw halves),
// the '+ denominator - 1' folds into the lea that adds them, so the rounding
// is free. Computing it from the remainder instead -- quotient plus
// (remainder != 0), which is total and needs no precondition -- measured
// +2 to +4 % user CPU on the CLAUDE.md 'd = 2' benchmark, because that
// version has to materialise the remainder with a multiply and a subtract.
// The precondition above is the price of the free version.
//
// (No power-of-two special case here, and vsearch's round_up_to_8() is not
//  one either way -- it lives in that project, not this one. For a constant
//  power-of-two denominator GCC already emits the bitmask by itself:
//  ((x + 3) / 4) * 4 and (x + 3) & ~3 compile to the same two instructions,
//  checked at -O3. Writing it by hand would only pay for a run-time
//  denominator known to be a power of two, and swarm has none.)


// C++14 refactoring: constexpr
// C++17 refactoring: [[nodiscard]]
// C++20 refactoring: std::unsigned_integral concept, in place of the
// type parameter and the static_assert below
template <typename Unsigned = std::uint64_t>
auto ceil_divide(Unsigned const numerator,
                 Unsigned const denominator) noexcept -> Unsigned {
  static_assert(std::is_unsigned<Unsigned>::value,
                "ceil_divide requires an unsigned integer type");
  assert(denominator >= 1);
  assert(numerator <= std::numeric_limits<Unsigned>::max() - (denominator - 1));
  return (numerator + denominator - 1) / denominator;
}


// refactoring: C++14 tests

/*

// test return type:
static_assert(std::is_same<decltype(ceil_divide<std::uint64_t>(0, 1)),
              std::uint64_t>::value, "");

// test default type parameter:
static_assert(ceil_divide(0ULL, 1ULL) == ceil_divide<std::uint64_t>(0, 1), "");

// test return values
static_assert(ceil_divide<std::uint64_t>(0, 1) == 0, "");
static_assert(ceil_divide<std::uint64_t>(0, 8) == 0, "");
static_assert(ceil_divide<std::uint64_t>(1, 8) == 1, "");
static_assert(ceil_divide<std::uint64_t>(7, 8) == 1, "");
static_assert(ceil_divide<std::uint64_t>(8, 8) == 1, "");
static_assert(ceil_divide<std::uint64_t>(9, 8) == 2, "");
static_assert(ceil_divide<std::uint64_t>(11, 8) == 2, "");
static_assert(ceil_divide<std::uint64_t>(16, 8) == 2, "");
static_assert(ceil_divide<std::uint64_t>(17, 8) == 3, "");

// the qgram_diff call site (denominator = 2 * qgramlength = 10)
static_assert(ceil_divide<std::uint64_t>(0, 10) == 0, "");
static_assert(ceil_divide<std::uint64_t>(1, 10) == 1, "");
static_assert(ceil_divide<std::uint64_t>(10, 10) == 1, "");
static_assert(ceil_divide<std::uint64_t>(11, 10) == 2, "");

*/

#endif
