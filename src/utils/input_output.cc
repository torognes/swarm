/*
    SWARM

    Copyright (C) 2012-2025 Torbjorn Rognes and Frederic Mahe

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

#include <cstdio>  // FILE, fdopen
#include <cstring>  // strcmp
#include <unistd.h>  // dup, STDIN_FILENO, STDOUT_FILENO


auto fopen_input(const char * filename) -> std::FILE *
{
  /* open the input stream given by filename, but use stdin if name is - */
  std::FILE * input_stream = nullptr;

  if (std::strcmp(filename, "-") == 0) {
    const int file_descriptor = dup(STDIN_FILENO);
    input_stream = file_descriptor > 0 ? fdopen(file_descriptor, "rb") : nullptr;
  }
  else {
    input_stream = fopen(filename, "rb");  // refactoring: prefer std::fstream
  }

  return input_stream;
}


auto fopen_output(const char * filename) -> std::FILE *
{
  /* open the output stream given by filename, but use stdout if name is - */
  std::FILE * output_stream {nullptr};

  if (std::strcmp(filename, "-") == 0) {
    const int file_descriptor = dup(STDOUT_FILENO);
    output_stream = file_descriptor > 0 ? fdopen(file_descriptor, "w") : nullptr;
  }
  else {
    output_stream = fopen(filename, "w");  // refactoring: prefer std::fstream
  }

  return output_stream;
}
