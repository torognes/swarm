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

#ifndef HASHER_GENERIC_H
#define HASHER_GENERIC_H

#include <cstdlib>  // std::size_t


// inspired by a talk given by Victor Ciura (2024) So You Think You
// Can Hash, CppCon 2024

// generic wrapper that turns any conforming hashing algorithm into a
// callable object compatible with standard hashing
template <class HashAlgorithm>
struct GenericHash {
  using result_type = typename HashAlgorithm::result_type;

  template <class Type>
  auto operator()(Type const & type) const noexcept -> result_type {
    HashAlgorithm hasher;
    hash_append(hasher, type);
    return static_cast<result_type>(hasher);
  }
};

// usage: std::unordered_set<View<char>, GenericHash<fnv1a>> headers;


// hash_append() overloads for primitive types: (unused for now)
//
// template <class HashAlgorithm>
// auto hash_append(HashAlgorithm & hasher, int integer) -> void {
//     hasher(&integer, sizeof(integer));
// }

// template <class HashAlgorithm, class Type>
// auto hash_append(HashAlgorithm & hasher, Type * ptr) -> void {
//   hasher(&ptr, sizeof(ptr));
// }

// hash_append() overloads for complex types:
// - see view.hpp


// tests

// auto main() -> int {
//     char const* str = "test";
//     auto const view  = View<char>{str, 4};
//     auto const view2 = View<char>{str, 3};

//     std::unordered_set<View<char>, GenericHash<fnv1a>> headers;
//     headers.insert(view);
//     std::cout << headers.size() << '\n';
//     headers.insert(view);
//     std::cout << headers.size() << '\n';
//     headers.insert(view2);
//     std::cout << headers.size() << '\n';

//     return 0;
// }


#endif // HASHER_GENERIC_H
