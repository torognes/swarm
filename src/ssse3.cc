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


#ifdef __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#endif

#ifdef __SSSE3__

/*
  SSSE3 specific code for x86-64

  Only include if __SSSE3__ is defined, which is done by the
  gcc compiler when the -mssse3 option or similar is given.

  This code requires the _mm_shuffle_epi8 intrinsic implemented
  with the PSHUFB instruction on the CPU. That instruction was
  available starting with the Intel Core architecture in 2006.
*/

#include <tmmintrin.h>  // _mm_shuffle_epi8

using WORD = unsigned short;
using BYTE = unsigned char;


/* 8-bit version with 16 channels */

auto dprofile_shuffle8(BYTE * dprofile,
                       BYTE * score_matrix,
                       BYTE * dseq_byte) -> void
{
  auto dprofile_m128i = reinterpret_cast<__m128i *>(dprofile);
  auto score_matrix_m128i = reinterpret_cast<__m128i *>(score_matrix);
  auto dseq_m128i = reinterpret_cast<__m128i *>(dseq_byte);
  const auto m0 = _mm_load_si128(dseq_m128i + 0);
  const auto m1 = _mm_load_si128(dseq_m128i + 1);
  const auto m2 = _mm_load_si128(dseq_m128i + 2);
  const auto m3 = _mm_load_si128(dseq_m128i + 3);
  __m128i t1;
  __m128i t2;
  __m128i t3;
  __m128i t4;

  auto profline8 = [&](const long long int j) {
    const auto t0 = _mm_load_si128(score_matrix_m128i + 2 * j);
    t1 = _mm_shuffle_epi8(t0, m0);
    t2 = _mm_shuffle_epi8(t0, m1);
    t3 = _mm_shuffle_epi8(t0, m2);
    t4 = _mm_shuffle_epi8(t0, m3);

    _mm_store_si128(dprofile_m128i + 4 * j + 0, t1);
    _mm_store_si128(dprofile_m128i + 4 * j + 1, t2);
    _mm_store_si128(dprofile_m128i + 4 * j + 2, t3);
    _mm_store_si128(dprofile_m128i + 4 * j + 3, t4);
  };

  profline8(0);  // -/gap/no nucleotide (0)
  profline8(1);  // A (1)
  profline8(2);  // C (2)
  profline8(3);  // G (3)
  profline8(4);  // T (4)
}


/* 16-bit version with 8 channels */

void dprofile_shuffle16(WORD * dprofile,
                        WORD * score_matrix,
                        BYTE * dseq_byte)
{
  static constexpr unsigned int channels {8};  // does 8 represent the number of channels?
  auto dprofile_m128i = reinterpret_cast<__m128i *>(dprofile);
  auto score_matrix_m128i = reinterpret_cast<__m128i *>(score_matrix);
  auto dseq_m128i = reinterpret_cast<__m128i *>(dseq_byte);

  const __m128i zero = _mm_setzero_si128();
  const __m128i one = _mm_set1_epi16(1);

  const __m128i t0 = _mm_load_si128(dseq_m128i + 0);

  __m128i m0 = _mm_unpacklo_epi8(t0, zero);
  m0 = _mm_slli_epi16(m0, 1);
  __m128i t1 = _mm_adds_epu16(m0, one);
  t1 = _mm_slli_epi16(t1, channels);
  m0 = _mm_or_si128(m0, t1);

  __m128i m1 = _mm_unpackhi_epi8(t0, zero);
  m1 = _mm_slli_epi16(m1, 1);
  __m128i t2 = _mm_adds_epu16(m1, one);
  t2 = _mm_slli_epi16(t2, channels);
  m1 = _mm_or_si128(m1, t2);

  const __m128i t3 = _mm_load_si128(dseq_m128i + 1);

  __m128i m2 = _mm_unpacklo_epi8(t3, zero);
  m2 = _mm_slli_epi16(m2, 1);
  __m128i t4 = _mm_adds_epu16(m2, one);
  t4 = _mm_slli_epi16(t4, channels);
  m2 = _mm_or_si128(m2, t4);

  __m128i m3 = _mm_unpackhi_epi8(t3, zero);
  m3 = _mm_slli_epi16(m3, 1);
  __m128i t5 = _mm_adds_epu16(m3, one);
  t5 = _mm_slli_epi16(t5, channels);
  m3 = _mm_or_si128(m3, t5);

  __m128i u0;
  __m128i u1;
  __m128i u2;
  __m128i u3;
  __m128i u4;

  auto profline16 = [&](const long long int j) {
    u0 = _mm_load_si128(score_matrix_m128i + 4 * j);  // refactoring: const?
    u1 = _mm_shuffle_epi8(u0, m0);
    u2 = _mm_shuffle_epi8(u0, m1);
    u3 = _mm_shuffle_epi8(u0, m2);
    u4 = _mm_shuffle_epi8(u0, m3);

    _mm_store_si128(dprofile_m128i + 4 * j + 0, u1);
    _mm_store_si128(dprofile_m128i + 4 * j + 1, u2);
    _mm_store_si128(dprofile_m128i + 4 * j + 2, u3);
    _mm_store_si128(dprofile_m128i + 4 * j + 3, u4);
  };

  profline16(0);
  profline16(1);
  profline16(2);
  profline16(3);
  profline16(4);
}

#else
#error __SSSE3__ not defined
#endif
#endif
