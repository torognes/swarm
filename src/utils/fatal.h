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

#include <iostream>


/* bare-exit form: terminate with a final newline (no error prefix) */
auto fatal() -> void;


/* message and exit with an error (variadic template with compile-time recursion) */

// refactoring C++11: constrain with type traits?
// refactoring C++17: use fold expression
// refactoring C++20: use concept "Printable"

namespace fatal_detail {

    // recursion base case: defined out-of-line in fatal.cc to avoid
    // multiple-definitions at link time
    auto print_then_exit() -> void;

    // recursive case: consume arguments one-by-one
    template<typename T, typename... Tail>
    auto print_then_exit(T head, Tail... tail) -> void {
        std::cerr << head;
        print_then_exit(tail...);
    }

}  // namespace fatal_detail


// public variadic: auto-prefix "\nError: " then forward
template<typename T, typename... Tail>
auto fatal(T head, Tail... tail) -> void {
    std::cerr << "\nError: ";
    fatal_detail::print_then_exit(head, tail...);
}
