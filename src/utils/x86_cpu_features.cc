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

#include "../swarm.h"
#include "fatal.h"
#include <cstdio>  // fprintf


#ifdef __x86_64__

#include <cpuid.h>  // __get_cpuid, __get_cpuid_count, bit_* masks

auto cpu_features_detect(struct Parameters & parameters) -> void
{
  // CPU registers:
  unsigned int eax {0};
  unsigned int ebx {0};
  unsigned int ecx {0};
  unsigned int edx {0};

  // leaf 1: standard feature flags
  if (__get_cpuid(1, &eax, &ebx, &ecx, &edx) == 0) {
    return;
  }
  parameters.mmx_present    = ((edx & bit_MMX)    != 0U) ? 1 : 0;
  parameters.sse_present    = ((edx & bit_SSE)    != 0U) ? 1 : 0;
  parameters.sse2_present   = ((edx & bit_SSE2)   != 0U) ? 1 : 0;
  parameters.sse3_present   = ((ecx & bit_SSE3)   != 0U) ? 1 : 0;
  parameters.ssse3_present  = ((ecx & bit_SSSE3)  != 0U) ? 1 : 0;
  parameters.sse41_present  = ((ecx & bit_SSE4_1) != 0U) ? 1 : 0;
  parameters.sse42_present  = ((ecx & bit_SSE4_2) != 0U) ? 1 : 0;
  parameters.popcnt_present = ((ecx & bit_POPCNT) != 0U) ? 1 : 0;
  parameters.avx_present    = ((ecx & bit_AVX)    != 0U) ? 1 : 0;

  // leaf 7, sub-leaf 0: extended feature flags
  static constexpr unsigned int extended_features_leaf {7};
  static constexpr unsigned int extended_features_subleaf {0};
  if (__get_cpuid_count(extended_features_leaf, extended_features_subleaf,
                        &eax, &ebx, &ecx, &edx) == 0) {
    return;
  }
  parameters.avx2_present   = ((ebx & bit_AVX2)   != 0U) ? 1 : 0;
}

auto cpu_features_test(struct Parameters & parameters) -> void {
  if (parameters.sse2_present == 0) {
    fatal(error_prefix, "This program requires a processor with SSE2 instructions.");
  }

  if (parameters.opt_disable_sse3)
    {
      parameters.sse3_present = 0;
      parameters.ssse3_present = 0;
      parameters.sse41_present = 0;
      parameters.sse42_present = 0;
      parameters.popcnt_present = 0;
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
