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

#include "input_output.hpp"
#include "fatal.hpp"  // fatal
#include <cstdio>  // fopen, FILE, fdopen, fclose, ferror
#include <string>  // std::string
#include <unistd.h>  // dup, STDIN_FILENO, STDOUT_FILENO


auto CloseFileHandle::operator()(std::FILE * file_handle) const -> void {
  // Two separate questions, both asked: std::ferror reports any read or
  // write that failed on this stream at any point, because stdio latches
  // the flag; std::fclose reports the final buffer flush, which can fail on
  // its own after every individual write appeared to succeed.
  bool const stream_failed = std::ferror(file_handle) != 0;
  bool const close_failed = std::fclose(file_handle) != 0;
  if (stream_failed or close_failed) {
    fatal("I/O error on a swarm file; the output may be incomplete.");
  }
}


namespace {

  // dup() returns -1 on failure; any non-negative value (including 0,
  // when a standard descriptor has been closed) is a valid descriptor.
  constexpr int invalid_fd {-1};

  auto is_dash(std::string const & filename) -> bool {
    return filename == "-";
  }

}  // end of anonymous namespace


auto fopen_input(std::string const & filename) -> FileHandle {
  /* open the input stream given by filename, but use stdin if name is - */
  std::FILE * input_stream {nullptr};

  if (is_dash(filename)) {
    auto const file_descriptor = dup(STDIN_FILENO);
    input_stream = file_descriptor != invalid_fd ? fdopen(file_descriptor, "rb") : nullptr;
  }
  else {
    input_stream = std::fopen(filename.c_str(), "rb");
  }

  return FileHandle{input_stream};
}


auto fopen_output(std::string const & filename) -> FileHandle {
  /* open the output stream given by filename, but use stdout if name is - */
  std::FILE * output_stream {nullptr};

  if (is_dash(filename)) {
    auto const file_descriptor = dup(STDOUT_FILENO);
    output_stream = file_descriptor != invalid_fd ? fdopen(file_descriptor, "w") : nullptr;
  }
  else {
    output_stream = std::fopen(filename.c_str(), "w");
  }

  return FileHandle{output_stream};
}
