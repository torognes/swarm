/*
    SWARM

    Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe

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

/* Please note: This code requires the pshufb instruction, which is
   part of the SSSE3 instruction set for x86_64 cpus. It seems not to be
   present on AMD cpus and not on the first generation Intel x86_64 cpus. */

/* 8-bit version with 16 channels */

void dprofile_shuffle8(BYTE * dprofile,
                       BYTE * score_matrix,
                       BYTE * dseq_byte)
{
  __m128i m0, m1, m2, m3, t0, t1, t2, t3, t4;
  __m128i * dseq = (__m128i*) dseq_byte;

  m0 = _mm_load_si128(dseq);
  m1 = _mm_load_si128(dseq+1);
  m2 = _mm_load_si128(dseq+2);
  m3 = _mm_load_si128(dseq+3);

#define profline8(j)                                    \
  t0 = _mm_load_si128((__m128i*)(score_matrix)+2*j);    \
  t1 = _mm_shuffle_epi8(t0, m0);                        \
  t2 = _mm_shuffle_epi8(t0, m1);                        \
  t3 = _mm_shuffle_epi8(t0, m2);                        \
  t4 = _mm_shuffle_epi8(t0, m3);                        \
  _mm_store_si128((__m128i*)(dprofile)+4*j+0, t1);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+1, t2);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+2, t3);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+3, t4)

  profline8(0);
  profline8(1);
  profline8(2);
  profline8(3);
  profline8(4);
}


/* 16-bit version with 8 channels */

void dprofile_shuffle16(WORD * dprofile,
                        WORD * score_matrix,
                        BYTE * dseq_byte)
{
  __m128i m0, m1, m2, m3;
  __m128i t0, t1, t2, t3, t4, t5;
  __m128i u0, u1, u2, u3, u4;
  __m128i zero, one;
  __m128i * dseq = (__m128i*) dseq_byte;

  zero = _mm_setzero_si128();
  one  = _mm_set1_epi16(1);

  t0 = _mm_load_si128(dseq+0);

  m0 = _mm_unpacklo_epi8(t0, zero);
  m0 = _mm_slli_epi16(m0, 1);
  t1 = _mm_adds_epu16(m0, one);
  t1 = _mm_slli_epi16(t1, 8);
  m0 = _mm_or_si128(m0, t1);

  m1 = _mm_unpackhi_epi8(t0, zero);
  m1 = _mm_slli_epi16(m1, 1);
  t2 = _mm_adds_epu16(m1, one);
  t2 = _mm_slli_epi16(t2, 8);
  m1 = _mm_or_si128(m1, t2);

  t3 = _mm_load_si128(dseq+1);

  m2 = _mm_unpacklo_epi8(t3, zero);
  m2 = _mm_slli_epi16(m2, 1);
  t4 = _mm_adds_epu16(m2, one);
  t4 = _mm_slli_epi16(t4, 8);
  m2 = _mm_or_si128(m2, t4);

  m3 = _mm_unpackhi_epi8(t3, zero);
  m3 = _mm_slli_epi16(m3, 1);
  t5 = _mm_adds_epu16(m3, one);
  t5 = _mm_slli_epi16(t5, 8);
  m3 = _mm_or_si128(m3, t5);

#define profline16(j)                                   \
  u0 = _mm_load_si128((__m128i*)(score_matrix)+4*j);    \
  u1 = _mm_shuffle_epi8(u0, m0);                        \
  u2 = _mm_shuffle_epi8(u0, m1);                        \
  u3 = _mm_shuffle_epi8(u0, m2);                        \
  u4 = _mm_shuffle_epi8(u0, m3);                        \
  _mm_store_si128((__m128i*)(dprofile)+4*j+0, u1);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+1, u2);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+2, u3);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+3, u4)

  profline16(0);
  profline16(1);
  profline16(2);
  profline16(3);
  profline16(4);
}

