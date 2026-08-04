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

#include "cpu_features.hpp"
#include "../../swarm.hpp"
#include "../../utils/fatal.hpp"
#include <cpuid.h>  // __cpuid_count, bit_* feature masks
#include "../../utils/print_view.hpp"  // fprint
#include <array>
#include <cstdint>  // int64_t
#include <cstdio>  // fputs()

namespace {
// Groups the four output registers of a CPUID query so each leaf can be
// returned as one named, const result instead of reusing shared variables.
struct cpuid_registers {
  unsigned int eax {0};
  unsigned int ebx {0};
  unsigned int ecx {0};
  unsigned int edx {0};
};

// All call sites query sub-leaf 0, so the sub-leaf is fixed here rather
// than passed in (avoids two adjacent same-type parameters). __cpuid_count
// has shipped in <cpuid.h> since GCC 4.4.
auto get_cpuid(unsigned int const leaf) noexcept -> cpuid_registers {
  cpuid_registers registers {};
  __cpuid_count(leaf, 0U, registers.eax, registers.ebx, registers.ecx, registers.edx);
  return registers;
}

// Read the low 32 bits of XCR0 via XGETBV. Must only be called when
// CPUID reports OSXSAVE, otherwise XGETBV raises #UD. cpu_features.cpp is
// compiled with the baseline target (no -mxsave), so the _xgetbv
// intrinsic is unavailable; use the equivalent one-instruction asm.
auto read_xcr0() noexcept -> unsigned int {
  unsigned int xcr0_lo {0};
  unsigned int xcr0_hi {0};
  __asm__ __volatile__("xgetbv" : "=a"(xcr0_lo), "=d"(xcr0_hi) : "c"(0U));
  static_cast<void>(xcr0_hi);
  return xcr0_lo;
}
}  // namespace

auto cpu_features_detect(struct Parameters & parameters) -> void
{
  // Feature masks (bit_MMX, bit_SSE, ...) come from <cpuid.h>. bit_OSXSAVE
  // is not defined by older <cpuid.h> versions (GCC 4.x), so spell it out.
  static constexpr unsigned int basic_leaf_mask {0xffU};    // CPUID.0:EAX low byte
  static constexpr unsigned int extended_features_leaf {7U};
  static constexpr unsigned int bit_osxsave {0x08000000U};  // CPUID.1:ECX bit 27
  static constexpr unsigned int xcr0_avx_state {0x6U};      // XMM | YMM

  // leaf 0: highest supported standard leaf, used to gate the queries below
  cpuid_registers const leaf0 = get_cpuid(0U);
  unsigned int const max_level = leaf0.eax & basic_leaf_mask;
  if (max_level < 1U) {
    return;
  }

  // leaf 1: standard feature flags
  cpuid_registers const leaf1 = get_cpuid(1U);
  parameters.mmx_present    = ((leaf1.edx & bit_MMX)    != 0U) ? 1 : 0;
  parameters.sse_present    = ((leaf1.edx & bit_SSE)    != 0U) ? 1 : 0;
  parameters.sse2_present   = ((leaf1.edx & bit_SSE2)   != 0U) ? 1 : 0;
  parameters.sse3_present   = ((leaf1.ecx & bit_SSE3)   != 0U) ? 1 : 0;
  parameters.ssse3_present  = ((leaf1.ecx & bit_SSSE3)  != 0U) ? 1 : 0;
  parameters.sse41_present  = ((leaf1.ecx & bit_SSE4_1) != 0U) ? 1 : 0;
  parameters.sse42_present  = ((leaf1.ecx & bit_SSE4_2) != 0U) ? 1 : 0;
  parameters.popcnt_present = ((leaf1.ecx & bit_POPCNT) != 0U) ? 1 : 0;

  // AVX/AVX2 are only usable if the OS has enabled saving of the YMM
  // register state: CPUID must report OSXSAVE and XCR0 (read via XGETBV)
  // must have both the SSE (bit 1) and AVX (bit 2) state-enable bits set.
  // Without this check an AVX-capable CPU on an old OS would be
  // over-reported.
  bool const avx_os_enabled =
    ((leaf1.ecx & bit_osxsave) != 0U) and ((read_xcr0() & xcr0_avx_state) == xcr0_avx_state);
  parameters.avx_present =
    (((leaf1.ecx & bit_AVX) != 0U) and avx_os_enabled) ? 1 : 0;

  if (max_level < extended_features_leaf) {
    return;
  }

  // leaf 7, sub-leaf 0: extended feature flags
  cpuid_registers const leaf7 = get_cpuid(extended_features_leaf);
  parameters.avx2_present =
    (((leaf7.ebx & bit_AVX2) != 0U) and avx_os_enabled) ? 1 : 0;
}

auto cpu_features_test(struct Parameters & parameters) -> void {
  if (parameters.sse2_present == 0) {
    fatal("This program requires a processor with SSE2 instructions.");
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
  // One entry per reported feature, in the order they are printed. A table
  // rather than ten near-identical if blocks: the printed order becomes a
  // property of the data instead of of the control flow, and adding a
  // feature is one line. The member is reached through a pointer-to-member
  // so that the flag and its name cannot drift apart.
  struct Feature {
    int64_t Parameters::* flag;
    char const * name;
  };

  static constexpr std::array<Feature, 10> features {{
      {&Parameters::mmx_present,    " mmx"},
      {&Parameters::sse_present,    " sse"},
      {&Parameters::sse2_present,   " sse2"},
      {&Parameters::sse3_present,   " sse3"},
      // Supplemental SSE3, introduced in 2006
      {&Parameters::ssse3_present,  " ssse3"},
      {&Parameters::sse41_present,  " sse4.1"},
      {&Parameters::sse42_present,  " sse4.2"},
      {&Parameters::popcnt_present, " popcnt"},
      {&Parameters::avx_present,    " avx"},
      {&Parameters::avx2_present,   " avx2"},
    }};

  fprint(parameters.logfile, "CPU features:     ");
  for (auto const & feature : features) {
    if ((parameters.*(feature.flag)) != 0) {
      static_cast<void>(std::fputs(feature.name, parameters.logfile));
    }
  }
  fprint(parameters.logfile, '\n');
}
