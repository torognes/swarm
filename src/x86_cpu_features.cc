/*
    SWARM

    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

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

// set to null if not x86-64
int64_t mmx_present {0};
int64_t sse_present {0};
int64_t sse2_present {0};
int64_t sse3_present {0};
int64_t ssse3_present {0};
int64_t sse41_present {0};
int64_t sse42_present {0};
int64_t popcnt_present {0};
int64_t avx_present {0};
int64_t avx2_present {0};

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

void cpu_features_detect()
{
  constexpr unsigned int post_pentium {7};  // new cpus: a & 0xff > 6
  constexpr unsigned int bit_mmx {23};
  constexpr unsigned int bit_sse {25};
  constexpr unsigned int bit_sse2 {26};
  constexpr unsigned int bit_sse3 {0};
  constexpr unsigned int bit_ssse3 {9};
  constexpr unsigned int bit_sse41 {19};
  constexpr unsigned int bit_sse42 {20};
  constexpr unsigned int bit_popcnt {23};
  constexpr unsigned int bit_avx {28};
  constexpr unsigned int bit_avx2 {5};

  unsigned int a {0};
  unsigned int b {0};
  unsigned int c {0};
  unsigned int d {0};

  cpuid(0, 0, a, b, c, d);
  unsigned int maxlevel = a & UINT8_MAX;

  if (maxlevel >= 1)
  {
    cpuid(1, 0, a, b, c, d);
    mmx_present    = (d >> bit_mmx) & 1;
    sse_present    = (d >> bit_sse) & 1;
    sse2_present   = (d >> bit_sse2) & 1;
    sse3_present   = (c >> bit_sse3) & 1;
    ssse3_present  = (c >> bit_ssse3) & 1;
    sse41_present  = (c >> bit_sse41) & 1;
    sse42_present  = (c >> bit_sse42) & 1;
    popcnt_present = (c >> bit_popcnt) & 1;
    avx_present    = (c >> bit_avx) & 1;

    if (maxlevel >= post_pentium)
    {
      cpuid(post_pentium, 0, a, b, c, d);
      avx2_present   = (b >> bit_avx2) & 1;
    }
  }
}

void cpu_features_show()
{
  fprintf(logfile, "CPU features:     ");
  if (mmx_present != 0){
    fprintf(logfile, " mmx");
  }
  if (sse_present != 0) {
    fprintf(logfile, " sse");
  }
  if (sse2_present != 0) {
    fprintf(logfile, " sse2");
  }
  if (sse3_present != 0) {
    fprintf(logfile, " sse3");
  }
  if (ssse3_present != 0) {
    fprintf(logfile, " ssse3"); // Supplemental SSSE3, introduced in 2006
  }
  if (sse41_present != 0) {
    fprintf(logfile, " sse4.1");
  }
  if (sse42_present != 0) {
    fprintf(logfile, " sse4.2");
  }
  if (popcnt_present != 0) {
    fprintf(logfile, " popcnt");
  }
  if (avx_present != 0) {
    fprintf(logfile, " avx");
  }
  if (avx2_present != 0) {
    fprintf(logfile, " avx2");
  }
  fprintf(logfile, "\n");
}

#endif
