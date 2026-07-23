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

#include <cstddef>  // std::size_t
#include <iostream>
#include <utility>


/* bare-exit form: terminate with a final newline (no error prefix) */
[[noreturn]] auto fatal() -> void;


/* message and exit with an error (variadic template with compile-time recursion) */

// refactoring C++11: constrain with type traits?
// refactoring C++17: use fold expression
// refactoring C++20: use concept "Printable"

namespace fatal_detail {

    // recursion base case: defined out-of-line in fatal.cpp to avoid
    // multiple-definitions at link time
    [[noreturn]] auto print_then_exit() -> void;

    // Explicit array-to-pointer conversion used right before streaming.
    // Without it, every fatal("...", ...) call site triggers the
    // cppcoreguidelines-pro-bounds-array-to-pointer-decay check at the
    // ostream operator<<, since string literals reach this template as
    // 'char const (&)[N]'.
    template<typename Type>
    auto explicit_decay(Type && value) noexcept -> Type && {
        return std::forward<Type>(value);
    }
    // The C-array parameter is intentional: string literals reach this
    // overload as 'char const (&)[N]', and the whole point of the helper
    // is to turn that into a pointer explicitly. std::array can't match
    // a string literal, so the C-array parameter can't be avoided.
    template<typename Type, std::size_t Size>
    // NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
    auto explicit_decay(Type const (&array)[Size]) noexcept -> Type const * {
        return static_cast<Type const *>(array);
    }

    // recursive case: consume arguments one-by-one (forwarding
    // references avoid per-level copies of std::string and friends)
    template<typename Head, typename... Tail>
    [[noreturn]] auto print_then_exit(Head && head, Tail &&... tail) -> void {
        std::cerr << explicit_decay(std::forward<Head>(head));
        print_then_exit(std::forward<Tail>(tail)...);
    }

}  // namespace fatal_detail


// public variadic: auto-prefix "\nError: " then forward
template<typename Head, typename... Tail>
[[noreturn]] auto fatal(Head && head, Tail &&... tail) -> void {
    std::cerr << "\nError: ";
    fatal_detail::print_then_exit(std::forward<Head>(head),
                                  std::forward<Tail>(tail)...);
}
