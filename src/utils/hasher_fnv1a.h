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

#ifndef HASHER_FNV1A_HPP
#define HASHER_FNV1A_HPP


#include "view.h"
#include <cstdlib>  // std::size_t


// Fowler-Noll-Vo (FNV-1A 64-bit) hash function
//
class fnv1a {
private:
  // initialize internal state
  static constexpr auto FNV_offset_basis = std::size_t{14695981039346656037U};
  static constexpr auto FNV_prime = std::size_t{1099511628211U};
  std::size_t hash = FNV_offset_basis;

public:
  // consume bytes and update internal state
  auto operator()(void const * key, std::size_t const length) noexcept -> void {
    auto const * ptr = static_cast<unsigned char const *>(key);
    auto const bytes = View<unsigned char>{ptr, length};
    for (auto const & byte : bytes) {
      hash = (hash ^ byte) * FNV_prime;
    }
  }

  // finalize internal state to size_t (conversion operator)
  using result_type = decltype(hash);
  explicit operator result_type() const noexcept {
    return hash;
  }
};


// Notes:
// 1) it is possible to use std algorithm and a lambda to extract the actual hashing
// from the rest of the loop (it yields the same assembler code):
//
// #include <numeric>  // std::accumulate
// auto operation = [](std::size_t & accumulator, unsigned char const & byte) -> decltype(hash) {
//   return (accumulator ^ byte) * FNV_prime;
//  };
// hash = std::accumulate(bytes.cbegin(), bytes.cend(), hash, operation);


// tests:

// auto main() -> int {
// {
//     fnv1a hasher;
//     hasher("test", 4);
//     std::cout << static_cast<std::size_t>(hasher) << '\n';
//     assert(static_cast<std::size_t>(hasher) == 18007334074686647077ULL);
// }

// {
//     fnv1a hasher;
//     hasher("", 0);
//     std::cout << static_cast<std::size_t>(hasher) << '\n';
//     assert(static_cast<std::size_t>(hasher) == 14695981039346656037ULL);
// }

// {
//     fnv1a hasher;
//     hasher("Amazon Redshift", 15);
//     std::cout << static_cast<std::size_t>(hasher) << '\n';
//     assert(static_cast<std::size_t>(hasher) == 7783490368944507294ULL);
// }

// {
//     fnv1a hasher;
//     hasher("Hello, world!", 13);
//     std::cout << static_cast<std::size_t>(hasher) << '\n';
//     assert(static_cast<std::size_t>(hasher) == 4094109891673226228ULL);
// }

//     return 0;
// }

#endif // HASHER_FNV1A_HPP
