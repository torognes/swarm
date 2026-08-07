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


#include "element_order.hpp"  // element_order, element_less
#include <algorithm>  // std::equal, std::lexicographical_compare, std::min
#include <cassert>
#include <cstddef>  // std::ptrdiff_t
#include <cstdlib>  // std::size_t
#include <iterator> // std::prev, std::next
#include <type_traits>  // std::is_arithmetic, std::remove_cv

#ifndef NDEBUG
#include <limits>
#endif


// The hash_append friend below is the integration point with the
// HashAlgorithm protocol used by utils/hasher_generic.hpp
// (Howard-Hinnant-style universal hashing). To use a View<char> as
// a key in an std::unordered_set, the call site picks an algorithm
// explicitly, e.g.:
//   std::unordered_set<View<char>, GenericHash<fnv1a>> seen;


// const-only, non-owning view over a contiguous sequence of elements
// (vectors or arrays) of any type Type (except std::vector<bool>).
//
// inspired by std::span (C++20) but read-only: the underlying data
// cannot be modified through the view.

template <typename Type = char>
class View {
public:
  // Default-constructed view is empty (data() == nullptr, size() == 0).
  // Provided so that aggregates containing a View can themselves be
  // default-constructed (e.g. seqinfo_s held in std::vector::resize()).
  View() noexcept = default;

  explicit View(Type const * const start, std::size_t const length) noexcept
    : start_ {start},
      length_ {length} {
    assert((start != nullptr) or (length == 0));
    assert(length <= max_ptrdiff);
  }

  // Operators
  //
  // The three comparison members below are restricted to arithmetic element
  // types, and are noexcept because of it: for an arithmetic Type the ordering
  // bottoms out in a built-in comparison (see element_order.hpp), which cannot
  // throw, whereas an arbitrary Type's operator< can. The restriction is what
  // makes the promise honest.
  //
  // Deliberately checked per-member rather than at class scope: a member
  // function of a class template is instantiated only when used, so a View
  // over a non-arithmetic Type stays perfectly legal to declare, iterate and
  // index -- View<var_s> in algod1_network.cpp does exactly that -- and only
  // an attempt to *compare* such a view is an error.
  auto operator==(View<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a View requires an arithmetic element type");
    return size() == other.size()
      and std::equal(cbegin(), cend(), other.cbegin());
  }
  auto operator!=(View<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a View requires an arithmetic element type");
    return not (*this == other);
  }
  // Ordering goes through element_order (see element_order.hpp), so that a
  // View<char> orders its bytes as unsigned char, like std::strcmp and
  // std::string, rather than as a possibly-signed char.
  auto operator<(View<Type> const & other) const noexcept -> bool {
    static_assert(comparable, "comparing a View requires an arithmetic element type");
    return std::lexicographical_compare(cbegin(), cend(),
                                        other.cbegin(), other.cend(),
                                        element_less<Type>{});
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
  // Predicate behind the comparison members' static_assert above. remove_cv is
  // needed because std::is_arithmetic<char const> is false, and View<Type const>
  // is an ordinary read-only instantiation that must stay comparable.
  static constexpr bool comparable =
    std::is_arithmetic<typename std::remove_cv<Type>::type>::value;

#ifndef NDEBUG
  // Upper bounds for the debug-build assertions above. Kept private so that
  // including this header does not export these names into the global
  // namespace (where they could shadow, or be shadowed by, an unrelated
  // max_size / max_ptrdiff elsewhere).
  // C++17 refactoring: [[maybe_unused]]
  static constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
  static constexpr auto max_size = std::numeric_limits<std::size_t>::max();
#endif

  Type const * start_ {};
  std::size_t  length_ {};
};


// A View over a whole container, with the element type deduced so that
// call sites neither spell it out nor reach for data():
//
//   make_view(cigar_string_)
//
// rather than
//
//   View<char>{cigar_string_.data(), cigar_string_.size()}
//
// A prefix or an interior slice comes from composing with the members
// that already exist -- make_view(vec).first(count), or
// make_view(vec).subview(offset, count), which also carries the bounds
// assertions that an open-coded &vec[offset] cannot. That composition is
// why there is no offset/count overload here.
//
// noexcept: the containers this is called with (std::vector, std::array,
// std::string) all have noexcept data() and size(), so the only other
// operation left is View's noexcept constructor. A container whose
// accessors can throw is not a supported argument.
//
// constexpr: what C++11 requires of a constexpr function is the shape of
// the declaration, not that every call can be folded -- one return
// statement, a literal return type (View has an implicitly constexpr
// defaulted default constructor and two literal members), and literal
// parameter types (a reference type is one). No standard container has a
// constexpr data() before C++17, so under C++11 the specializations used
// here are simply not constant expressions; the keyword costs nothing and
// starts working the day the standard level moves.
template <typename Container>
constexpr auto make_view(Container const & container) noexcept
  -> View<typename Container::value_type> {
  return View<typename Container::value_type>{container.data(), container.size()};
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
