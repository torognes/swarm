/*
  SWARM

  Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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

#include <algorithm>
#include <cmath>
#include <type_traits>


// Computes the greatest common divisor of two integers
// refactoring C++17: replace with std::gcd()
template <typename Integral>
auto gcd(Integral num1, Integral num2) -> Integral {

  static_assert(std::is_integral<Integral>::value, "Integral required.");

  // gcd(a, b) = gcd(|a|, |b|)
  num1 = std::abs(num1);
  num2 = std::abs(num2);

  // gcd(a, b) = gcd(b, a) and gcd(a, 0) = gcd(0, a) = a
  // This covers the case gcd(0, 0) = 0
  if (num1 == 0) { return num2; }
  if (num2 == 0) { return num1; }

  while (num1 != 0) {
    if (num1 < num2) {
      std::swap(num1, num2);
    }
    num1 = num1 % num2;
  }
  return num2;
}

// refactoring C++20
// static_assert(gcd(0, 0) == 0);
// static_assert(gcd(0, 5) == 5);
// static_assert(gcd(5, 0) == 5);
// static_assert(gcd(5, 5) == 5);
// static_assert(gcd(-5, 0) == 5);
// static_assert(gcd(15, 9) == 3);
// static_assert(gcd(13, 6) == 1);   // one prime
// static_assert(gcd(6, 11) == 1);   // one prime
// static_assert(gcd(13, 11) == 1);  // two primes
// static_assert(gcd(100, 24) == 4);
