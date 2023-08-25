/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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

#include "swarm.h"
#include "utils/fatal.h"
#include <cstdint>  // int64_t


// set to null if not x86-64
int64_t ssse3_present {0};
int64_t sse41_present {0};
int64_t popcnt_present {0};

#ifdef __x86_64__

void cpuid(unsigned int f1,
           unsigned int f2,
           unsigned int & a,
           unsigned int & b,
           unsigned int & c,
           unsigned int & d);

void cpuid(unsigned int f1,
           unsigned int f2,
           unsigned int & a,
           unsigned int & b,
           unsigned int & c,
           unsigned int & d)
{
  __asm__ __volatile__ ("cpuid"
                        : "=a" (a), "=b" (b), "=c" (c), "=d" (d)
                        : "a" (f1), "c" (f2));
}

auto cpu_features_detect(struct Parameters & parameters) -> void
{
  static constexpr unsigned int post_pentium {7};  // new cpus: a & 0xff > 6
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

  unsigned int a {0};
  unsigned int b {0};
  unsigned int c {0};
  unsigned int d {0};

  cpuid(0, 0, a, b, c, d);
  const unsigned int maxlevel = a & UINT8_MAX;

  if (maxlevel >= 1)
  {
    cpuid(1, 0, a, b, c, d);
    parameters.mmx_present    = (d >> bit_mmx) & 1U;
    parameters.sse_present    = (d >> bit_sse) & 1U;
    parameters.sse2_present   = (d >> bit_sse2) & 1U;
    parameters.sse3_present   = (c >> bit_sse3) & 1U;
    ssse3_present  = (c >> bit_ssse3) & 1U;
    sse41_present  = (c >> bit_sse41) & 1U;
    parameters.sse42_present  = (c >> bit_sse42) & 1U;
    popcnt_present = (c >> bit_popcnt) & 1U;
    parameters.avx_present    = (c >> bit_avx) & 1U;

    if (maxlevel >= post_pentium)
    {
      cpuid(post_pentium, 0, a, b, c, d);
      parameters.avx2_present   = (b >> bit_avx2) & 1U;
    }
  }
}

auto cpu_features_test(struct Parameters & parameters) -> void {
  if (parameters.sse2_present == 0) {
    fatal(error_prefix, "This program requires a processor with SSE2 instructions.");
  }

  if (parameters.opt_disable_sse3)
    {
      parameters.sse3_present = 0;
      ssse3_present = 0;
      sse41_present = 0;
      parameters.sse42_present = 0;
      popcnt_present = 0;
      parameters.avx_present = 0;
      parameters.avx2_present = 0;
    }
}

auto cpu_features_show(struct Parameters const & parameters) -> void
{
  fprintf(logfile, "CPU features:     ");
  if (parameters.mmx_present != 0){
    fprintf(logfile, " mmx");
  }
  if (parameters.sse_present != 0) {
    fprintf(logfile, " sse");
  }
  if (parameters.sse2_present != 0) {
    fprintf(logfile, " sse2");
  }
  if (parameters.sse3_present != 0) {
    fprintf(logfile, " sse3");
  }
  if (ssse3_present != 0) {
    fprintf(logfile, " ssse3"); // Supplemental SSE3, introduced in 2006
  }
  if (sse41_present != 0) {
    fprintf(logfile, " sse4.1");
  }
  if (parameters.sse42_present != 0) {
    fprintf(logfile, " sse4.2");
  }
  if (popcnt_present != 0) {
    fprintf(logfile, " popcnt");
  }
  if (parameters.avx_present != 0) {
    fprintf(logfile, " avx");
  }
  if (parameters.avx2_present != 0) {
    fprintf(logfile, " avx2");
  }
  fprintf(logfile, "\n");
}

#endif
