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

#ifndef SWARM_UTILS_SPAN_H
#define SWARM_UTILS_SPAN_H


#include "view.hpp"
#include <algorithm>  // std::min
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <iterator> // std::prev, std::next

#ifndef NDEBUG
#include <limits>
#endif


// non-owning, mutable view over a contiguous sequence of elements
// (vectors or arrays) of any type Type (except std::vector<bool>).
//
// inspired by std::span (C++20). Mirrors View<Type> for the mutable
// case: elements may be modified through the span. Implicitly
// converts to View<Type> wherever read-only access is sufficient.

template <typename Type = char>
class Span {
public:
  // Default-constructed span is empty (data() == nullptr, size() == 0).
  Span() noexcept = default;

  explicit Span(Type * const start, std::size_t const length) noexcept
    : start_ {start},
      length_ {length} {
    assert((start != nullptr) or (length == 0));
    assert(length <= max_ptrdiff);
  }

  // Explicit conversion to read-only View<Type>: callers wanting to
  // hand a Span to a View-consuming API must opt in via
  // static_cast<View<Type>>(span) or View<Type>{span}.
  explicit operator View<Type>() const noexcept {
    return View<Type>{start_, length_};
  }

  // Iterators
  auto begin()  const noexcept -> Type * { return data(); }
  auto end() const noexcept -> Type * {
    auto const distance = static_cast<std::ptrdiff_t>(size());
    return std::next(data(), distance);
  }
  auto cbegin() const noexcept -> Type const * { return data(); }
  auto cend() const noexcept -> Type const * {
    auto const distance = static_cast<std::ptrdiff_t>(size());
    return std::next(static_cast<Type const *>(data()), distance);
  }
  auto rbegin() const noexcept -> std::reverse_iterator<Type *> {
    return std::reverse_iterator<Type *>(end());
  }
  auto rend() const noexcept -> std::reverse_iterator<Type *> {
    return std::reverse_iterator<Type *>(begin());
  }

  // Element access
  auto front() const noexcept -> Type & {
    assert(not empty());
    return *data();
  }
  auto back() const noexcept -> Type & {
    assert(not empty());
    return *std::prev(end());
  }
  auto data() const noexcept -> Type * { return start_; }
  auto operator[](std::size_t const index) const noexcept -> Type & {
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

  // Subspans
  auto subspan(std::size_t const offset, std::size_t const count) const noexcept -> Span {
    assert(offset <= size());
    assert(count <= size() - offset);
    auto const distance = static_cast<std::ptrdiff_t>(offset);
    auto * const new_start = std::next(data(), distance);
    return Span{new_start, count};
  }
  auto first(std::size_t const count) const noexcept -> Span {
    return subspan(0, count);
  }
  auto last(std::size_t const count) const noexcept -> Span {
    assert(count <= size());
    return subspan(size() - count, count);
  }
  auto drop(std::size_t const count) const noexcept -> Span {
    // drop n first items, return empty if n is >= size()
    auto const offset = std::min(count, size());
    assert(offset <= size());
    return subspan(offset, size() - offset);
  }

private:
#ifndef NDEBUG
  // C++17 refactoring: [[maybe_unused]]
  static constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
  static constexpr auto max_size = std::numeric_limits<std::size_t>::max();
#endif

  Type * start_ {};
  std::size_t length_ {};
};


#endif // SWARM_UTILS_SPAN_H
