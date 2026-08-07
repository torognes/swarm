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

// Portable XOR + popcount fallback. Opt-in via `make SWARM_PORTABLE_FALLBACK=1`,
// which swaps this file in for the ISA-specific arch/<isa>/qgram_compare.cpp.
// Intended as a building block for porting swarm to architectures that
// don't yet have a hand-written SIMD path.
//
// Note: this only replaces compareqgramvectors. The search8 / search16
// intrinsics under arch/<isa>/ have no portable equivalent yet, so a full
// port still requires per-ISA work for those.
//
// C++20 refactoring: __builtin_popcountll becomes std::popcount; the
// memcpy guards against strict-aliasing UB when the compiler is paranoid
// (modern GCC/Clang elide it to a single load on every supported target).

#include "../../utils/qgram_compare.hpp"
#include "../../utils/cpu_features.hpp"  // Cpu_features
#include "../../utils/qgram_array.hpp"   // qgramvectorbytes
#include <cstdint>  // uint64_t
#include <cstring>  // std::memcpy


auto compareqgramvectors(Qgram_vector const & lhs, Qgram_vector const & rhs,
                         Cpu_features const & cpu_features) -> uint64_t
{
  static_cast<void>(cpu_features);  // unused: portable path has no runtime dispatch
  static constexpr auto n_words = qgramvectorbytes / sizeof(uint64_t);  // 16
  uint64_t count {0};
  uint64_t lhs_word {0};
  uint64_t rhs_word {0};

  // data() rather than &lhs, unlike the SIMD kernels: this path copies the
  // bytes out instead of loading them as vectors, so it wants the buffer,
  // not an over-aligned pointer to reinterpret.
  auto const * lhs_bytes = lhs.data();
  auto const * rhs_bytes = rhs.data();

  // std::memcpy avoids the strict-aliasing undefined behaviour of punning
  // an unsigned char buffer through a uint64_t* (same reasoning as
  // variants.cpp nt_set); optimisers fold each one back into a single
  // load, so this is not a copy at run time.
  // C++20 refactoring: std::bit_cast
  for (auto i = 0ULL; i < n_words; ++i) {
    std::memcpy(&lhs_word, lhs_bytes + (i * sizeof(uint64_t)), sizeof(uint64_t));
    std::memcpy(&rhs_word, rhs_bytes + (i * sizeof(uint64_t)), sizeof(uint64_t));
    count += static_cast<uint64_t>(__builtin_popcountll(lhs_word ^ rhs_word));
  }

  return count;
}
