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

#ifndef SWARM_UTILS_THREAD_COUNT_H
#define SWARM_UTILS_THREAD_COUNT_H

#include <cassert>  // assert()
#include <cstddef>  // std::size_t
#include <cstdint>  // std::uint32_t


// The -t/--threads value, carrying the range it is allowed to hold: 1 to 512
// (man swarm, --threads). The command line reports out-of-range input with
// fatal() before building one of these, so the constructor asserts what the
// caller has already established rather than checking it again; the type
// exists to keep a validated count validated from then on.
//
// There is deliberately no conversion operator. count() is the only way out,
// and it returns std::size_t because that is what every consumer asks for:
// the ThreadRunner constructor and the per-thread vectors in Scanner and
// QgramDiffer are all sized by it. Without an implicit conversion a thread
// count cannot drift into arithmetic it was never meant to take part in.
class ThreadCount {
public:
  static constexpr std::uint32_t minimum {1};
  static constexpr std::uint32_t maximum {512};

  ThreadCount() = default;

  explicit ThreadCount(std::uint32_t const requested) noexcept
    : count_(requested) {
    assert(requested >= minimum);
    assert(requested <= maximum);
  }

  auto count() const noexcept -> std::size_t { return count_; }

private:
  std::uint32_t count_ {minimum};
};

#endif  // SWARM_UTILS_THREAD_COUNT_H
