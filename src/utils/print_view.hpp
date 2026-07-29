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

#ifndef SWARM_UTILS_PRINT_VIEW_H
#define SWARM_UTILS_PRINT_VIEW_H


#include "view.hpp"  // View<char>
#include <cstddef>  // std::size_t
#include <cstdio>  // std::FILE, std::fwrite


// Emit the bytes of a View to a stream.
//
// This is what a "%.*s" conversion was used for before: pass a counted run
// of characters that is not NUL-terminated (a slice of a header inside a
// shared buffer, a decoded sequence). std::fwrite takes both its element
// size and its count as an std::size_t, and the count is what View::size()
// already returns, so nothing is narrowed to the int that "%.*s" requires
// -- and there is no format string to parse at run time.
//
// Not exactly interchangeable with "%.*s": that conversion stops at an
// embedded NUL as well as at the precision, while fwrite always emits the
// full byte count. No header and no sequence reaching swarm's printers
// contains a NUL, because the reader truncates the line at the first one
// (Line_buffer::read_next), so the two agree on every input -- by
// construction of the input rather than by rejection of it. If swarm ever
// rejects NUL instead of truncating, this is the note to re-read.
//
// It lives here rather than in view.hpp so that View stays a pure data type
// and <cstdio> stays out of a header that nearly every translation unit
// reaches through seqinfo.hpp. Same reasoning as view_stream.hpp and
// <ostream>.

inline auto fprint(std::FILE * output_handle, View<char> const text) -> void
{
  // An empty view may carry a null pointer, and passing one to fwrite is
  // undefined even with a zero count. Reachable: a header whose abundance
  // annotation starts at offset 0 prints an empty leading slice.
  if (text.empty()) { return; }
  // fwrite's two size arguments are an element size and a count, in that
  // order, and transposing them compiles. Naming the element size says
  // which is which; sizeof(char) is 1 by definition, so this is about the
  // reader, not about a platform where it could differ.
  static constexpr std::size_t element_size = sizeof(char);
  static_cast<void>(std::fwrite(text.data(), element_size, text.size(), output_handle));
}

#endif // SWARM_UTILS_PRINT_VIEW_H
