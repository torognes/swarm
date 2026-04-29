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

#ifndef SWARM_UTILS_VIEW_H
#define SWARM_UTILS_VIEW_H


#include <algorithm>  // std::equal, std::lexicographical_compare, std::min
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdint>  // std::uint64_t
#include <cstdlib>  // std::size_t
#include <functional>  // std::hash
#include <iterator> // std::prev, std::next
#include <type_traits>  // std::is_arithmetic

#ifndef NDEBUG
#include <limits>
// C++17 refactoring: [[maybe_unused]]
constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
constexpr auto max_size = std::numeric_limits<std::size_t>::max();
#endif


// TODO:
//  - the hash_append friend below remains for use with custom hash
//    algorithms (Howard Hinnant style); std::hash<View<char>> at
//    the bottom of this header covers the std::unordered_set use
//    case.


// const-only, non-owning view over a contiguous sequence of elements
// (vectors or arrays) of any type Type (except std::vector<bool>).
//
// inspired by std::span (C++20) but read-only: the underlying data
// cannot be modified through the view.

template <typename Type = char>
class View {
public:
  explicit View(Type const * start, std::size_t const length) noexcept
    : start_ {start},
      length_ {length} {
    assert((start != nullptr) or (length == 0));
    assert(length <= max_ptrdiff);
  }

  // Operators
  auto operator==(View<Type> const & other) const noexcept -> bool {
    return size() == other.size()
      and std::equal(cbegin(), cend(), other.cbegin());
  }
  auto operator!=(View<Type> const & other) const noexcept -> bool {
    return not (*this == other);
  }
  auto operator<(View<Type> const & other) const noexcept -> bool {
    static_assert(std::is_arithmetic<Type>::value,
                  "View::operator< requires an arithmetic element type");
    return std::lexicographical_compare(cbegin(), cend(),
                                        other.cbegin(), other.cend());
  }

  // Iterators
  auto begin()  const noexcept -> Type const * { return data(); }
  auto cbegin() const noexcept -> Type const * { return data(); }
  auto end() const noexcept -> Type const * {
    auto const distance = static_cast<std::ptrdiff_t>(size());
    return std::next(data(), distance);
  }
  auto cend() const noexcept -> Type const * {
    return end();
  }
  auto rbegin() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(end());
  }
  auto crbegin() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(cend());
  }
  auto rend() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(begin());
  }
  auto crend() const noexcept -> std::reverse_iterator<Type const *> {
    return std::reverse_iterator<Type const *>(cbegin());
  }

  // Element access
  // C++17 refactoring: [[nodiscard]]
  auto front() const noexcept -> Type const & {
    assert(not empty());
    return *data();
  }
  auto back() const noexcept -> Type const & {
    assert(not empty());
    return *std::prev(end());
  }
  auto data() const noexcept -> Type const * { return start_; }
  auto operator[](std::size_t const index) const noexcept -> Type const & {
    assert(index < size());
    auto const distance = static_cast<std::ptrdiff_t>(index);
    return *std::next(data(), distance);
  }

  // Observers
  auto size() const noexcept -> std::size_t { return length_; }
  auto size_bytes() const noexcept -> std::size_t {
    assert(size() <= (max_size / sizeof(Type)));
    return size() * sizeof(Type);
  }
  auto empty() const noexcept -> bool { return size() == 0; }

  // Subviews
  auto subview(std::size_t const offset, std::size_t const count) const noexcept -> View {
    assert(offset <= size());
    assert(count <= size() - offset);
    auto const distance = static_cast<std::ptrdiff_t>(offset);
    auto const * new_start = std::next(data(), distance);
    return View{new_start, count};
  }
  auto first(std::size_t const count) const noexcept -> View {
    return subview(0, count);
  }
  auto last(std::size_t const count) const noexcept -> View {
    assert(count <= size());
    return subview(size() - count, count);
  }
  auto drop(std::size_t const count) const noexcept -> View {
    // drop n first items, return empty if n is >= size()
    auto const offset = std::min(count, size());
    assert(offset <= size());
    return subview(offset, size() - offset);
  }

  // hashing
  template <class HashAlgorithm>
  friend auto hash_append(HashAlgorithm & hasher, View<Type> const & view) noexcept -> void {
    hasher(view.data(), view.size_bytes());
  }

private:
  Type const * start_ {};
  std::size_t  length_ {};
};


// std::hash specialization, so that View<char> can be used as the
// key type of std::unordered_set / std::unordered_map.
//
// Implementation: FNV-1a 64-bit over the byte sequence. Adequate
// for short inputs such as fasta headers; not cryptographic. If
// adversarial input ever becomes a concern, swap for a stronger
// hash (e.g. SipHash with a runtime seed).
//
// Specialized for View<char> only on purpose: View<unsigned char>,
// View<std::uint64_t>, etc. would require byte reinterpretation
// machinery that is not yet justified by any caller.
namespace std {
  template <>
  struct hash<View<char>> {
    auto operator()(View<char> const & view) const noexcept -> std::size_t {
      static constexpr std::uint64_t fnv_offset_basis {14695981039346656037ULL};
      static constexpr std::uint64_t fnv_prime         {1099511628211ULL};
      std::uint64_t accumulator {fnv_offset_basis};
      for (auto const character : view) {
        accumulator ^= static_cast<std::uint64_t>(
                         static_cast<unsigned char>(character));
        accumulator *= fnv_prime;
      }
      return static_cast<std::size_t>(accumulator);
    }
  };
}


// tests:

// #include <algorithm>
// #include <vector>

// auto main() -> int {
//   std::vector<char> v = {'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'};
//   auto s = View<char>{v.data(), 5};
//   assert(s.size() == 5);
//   assert(s.size_bytes() == 5);
//   assert(! s.empty());
//   assert(s.front() == 'a');
//   assert(s.back() == 'e');
//   assert(s[1] == 'b');
//   for (auto c: s) {
//     printf("%c\n", c);
//   }
//   printf("\n");
//   auto report = [](char const &c) -> void { printf("%c\n", c); };
//   std::for_each(s.begin(), s.end(), report);
//   printf("\n");
//   // assert(s[5] == 'f');
//   auto s1 = View<char>{v.data(), 0};
//   std::for_each(s1.begin(), s1.end(), report);
//   printf("\n");
//   auto s2 = View<char>{v.data(), 10};
//   auto s3 = s2.first(2);
//   std::for_each(s3.begin(), s3.end(), report);
//   printf("\n");
//   auto s4 = s2.first(10);
//   std::for_each(s4.begin(), s4.end(), report);
//   printf("\n");
//   auto s5 = s2.first(0);
//   std::for_each(s5.begin(), s5.end(), report);
//   printf("\n");
//   auto s6 = s2.last(0);
//   std::for_each(s6.begin(), s6.end(), report);
//   printf("\n");
//   auto s7 = s2.last(2);
//   std::for_each(s7.begin(), s7.end(), report);
//   printf("\n");
//   std::for_each(s7.rbegin(), s7.rend(), report);
//   printf("\n");
//   std::for_each(s7.crbegin(), s7.crend(), report);
//   printf("\n");
// }

#endif // SWARM_UTILS_VIEW_H
