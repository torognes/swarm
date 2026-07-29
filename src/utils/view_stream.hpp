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

#ifndef SWARM_UTILS_VIEW_STREAM_H
#define SWARM_UTILS_VIEW_STREAM_H


#include "view.hpp"  // View<char>
#include <ios>  // std::streamsize
#include <ostream>  // std::ostream


// Insert the bytes of a View<char> into an output stream.
//
// This is what fatal() needs: it streams each of its arguments through
// std::cerr, so an argument has to be insertable. Without this operator a
// caller holding a view had to reach for the underlying pointer (and rely on
// a '\0' having been written after it) or materialize an std::string.
//
// std::ostream::write is used rather than a loop or an std::string, because
// the count is what View::size() already carries: the view need not be
// NUL-terminated, and no copy is made.
//
// It lives here rather than in view.hpp so that View stays a pure data type
// and <ostream> stays out of a header that nearly every translation unit
// reaches through seqinfo.hpp. Same reasoning as print_view.hpp and <cstdio>.

inline auto operator<<(std::ostream & output_stream, View<char> const & text) -> std::ostream &
{
  // An empty view may carry a null pointer; do not hand one to write(), even
  // with a zero count.
  if (text.empty()) { return output_stream; }
  output_stream.write(text.data(), static_cast<std::streamsize>(text.size()));
  return output_stream;
}

#endif // SWARM_UTILS_VIEW_STREAM_H
