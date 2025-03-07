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

#include "../swarm.h"
#include "fatal.h"
#include "x86_cpu_feature_popcnt.h"
#include "x86_cpu_feature_sse41.h"
#include "x86_cpu_feature_ssse3.h"
#include <cstdint>  // int64_t
#include <cstdio>  // fprintf
#include <limits>


// set to null if not x86-64
int64_t ssse3_present {0};
int64_t sse41_present {0};
int64_t popcnt_present {0};

#ifdef __x86_64__

// refactoring: rewrite using header 'cpuid.h' (GCC, clang)
auto cpuid(unsigned int leaf_level,
           unsigned int sublevel,
           unsigned int & eax,
           unsigned int & ebx,
           unsigned int & ecx,
           unsigned int & edx) -> void
{
  __asm__ __volatile__ ("cpuid"
                        : "=a" (eax), "=b" (ebx), "=c" (ecx), "=d" (edx)
                        : "a" (leaf_level), "c" (sublevel));
}

auto cpu_features_detect(struct Parameters & parameters) -> void
{
  static constexpr auto uint8_max = std::numeric_limits<uint8_t>::max();
  static constexpr unsigned int post_pentium {7};  // new cpus: eax & 0xff > 6
  static constexpr unsigned int bit_mmx {23};
  static constexpr unsigned int bit_sse {25};
  static constexpr unsigned int bit_sse2 {26};
  static constexpr unsigned int bit_sse3 {0};
  static constexpr unsigned int bit_ssse3 {9};
  static constexpr unsigned int bit_sse41 {19};
  static constexpr unsigned int bit_sse42 {20};
  static constexpr unsigned int bit_popcnt {23};
  static constexpr unsigned int bit_avx {28};
  static constexpr unsigned int bit_avx2 {5};

  // CPU registers:
  unsigned int eax {0};
  unsigned int ebx {0};
  unsigned int ecx {0};
  unsigned int edx {0};

  cpuid(0, 0, eax, ebx, ecx, edx);  // leaf 0
  const unsigned int maxlevel = eax & uint8_max;

  if (maxlevel == 0) {
    return;
  }

  cpuid(1, 0, eax, ebx, ecx, edx);  // leaf 1
  parameters.mmx_present    = (edx >> bit_mmx) & 1U;
  parameters.sse_present    = (edx >> bit_sse) & 1U;
  parameters.sse2_present   = (edx >> bit_sse2) & 1U;
  parameters.sse3_present   = (ecx >> bit_sse3) & 1U;
  parameters.ssse3_present = (ecx >> bit_ssse3) & 1U;
  ssse3_present = parameters.ssse3_present;
  parameters.sse41_present = (ecx >> bit_sse41) & 1U;
  sse41_present  = parameters.sse41_present;
  parameters.sse42_present  = (ecx >> bit_sse42) & 1U;
  parameters.popcnt_present = (ecx >> bit_popcnt) & 1U;
  popcnt_present = parameters.popcnt_present;
  parameters.avx_present    = (ecx >> bit_avx) & 1U;

  if (maxlevel >= post_pentium)
    {
      cpuid(post_pentium, 0, eax, ebx, ecx, edx);  // leaf 7
      parameters.avx2_present   = (ebx >> bit_avx2) & 1U;
    }
}

auto cpu_features_test(struct Parameters & parameters) -> void {
  if (parameters.sse2_present == 0) {
    fatal(error_prefix, "This program requires a processor with SSE2 instructions.");
  }

  if (parameters.opt_disable_sse3)
    {
      parameters.sse3_present = 0;
      parameters.ssse3_present = 0;
      ssse3_present = 0;
      parameters.sse41_present = 0;
      sse41_present = 0;
      parameters.sse42_present = 0;
      parameters.popcnt_present = 0;
      popcnt_present = 0;
      parameters.avx_present = 0;
      parameters.avx2_present = 0;
    }
}

auto cpu_features_show(struct Parameters const & parameters) -> void
{
  std::fprintf(parameters.logfile, "CPU features:     ");
  if (parameters.mmx_present != 0){
    std::fprintf(parameters.logfile, " mmx");
  }
  if (parameters.sse_present != 0) {
    std::fprintf(parameters.logfile, " sse");
  }
  if (parameters.sse2_present != 0) {
    std::fprintf(parameters.logfile, " sse2");
  }
  if (parameters.sse3_present != 0) {
    std::fprintf(parameters.logfile, " sse3");
  }
  if (parameters.ssse3_present != 0) {
    std::fprintf(parameters.logfile, " ssse3"); // Supplemental SSE3, introduced in 2006
  }
  if (parameters.sse41_present != 0) {
    std::fprintf(parameters.logfile, " sse4.1");
  }
  if (parameters.sse42_present != 0) {
    std::fprintf(parameters.logfile, " sse4.2");
  }
  if (parameters.popcnt_present != 0) {
    std::fprintf(parameters.logfile, " popcnt");
  }
  if (parameters.avx_present != 0) {
    std::fprintf(parameters.logfile, " avx");
  }
  if (parameters.avx2_present != 0) {
    std::fprintf(parameters.logfile, " avx2");
  }
  std::fprintf(parameters.logfile, "\n");
}

#endif
