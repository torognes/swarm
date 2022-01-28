/*
    SWARM

    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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

/* message and exit with an error (variadic template with compile-time recursion) */

// C++17: refactor with fold expression
// C++20: refactor with Printable concept

static char const * const error_prefix {"\nError: "};  // C++17: refactor by moving to fatal template

auto fatal() -> void;

template <typename T, typename... Tail>
auto fatal(T head, Tail... tail) -> void {
    std::cerr << head;
    fatal(tail...);
}
