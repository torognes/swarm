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
// refactoring? using VECTORTYPE = __m128i; (same as sse41?)


/* 8-bit version with 16 channels */

auto dprofile_shuffle8(BYTE * dprofile,
                       BYTE * score_matrix,
                       BYTE * dseq_byte) -> void
{
  auto * profile_db = reinterpret_cast<__m128i *>(dprofile);    // output
  auto * score_db = reinterpret_cast<__m128i *>(score_matrix);  // input
  auto * sequence_db = reinterpret_cast<__m128i *>(dseq_byte);  // input
  const auto seq_chunk0 = _mm_load_si128(sequence_db + 0);  // 16 nucleotides
  const auto seq_chunk1 = _mm_load_si128(sequence_db + 1);  // next 16
  const auto seq_chunk2 = _mm_load_si128(sequence_db + 2);  // next 16
  const auto seq_chunk3 = _mm_load_si128(sequence_db + 3);  // final 16 (total of 64)

  auto profline8 = [&](const long long int nuc) {
    // scores: 16 scores from the score matrix, matching the
    // nucleotide 'nuc'; five different nucleotides (0, 1, 2, 3, 4),
    // so five possible rows of scores
    const auto scores = _mm_load_si128(score_db + 2 * nuc);

    _mm_store_si128(profile_db + 4 * nuc + 0, _mm_shuffle_epi8(scores, seq_chunk0));
    _mm_store_si128(profile_db + 4 * nuc + 1, _mm_shuffle_epi8(scores, seq_chunk1));
    _mm_store_si128(profile_db + 4 * nuc + 2, _mm_shuffle_epi8(scores, seq_chunk2));
    _mm_store_si128(profile_db + 4 * nuc + 3, _mm_shuffle_epi8(scores, seq_chunk3));
  };

  profline8(0);  // -/gap/no nucleotide (0)
  profline8(1);  // A (1)
  profline8(2);  // C (2)
  profline8(3);  // G (3)
  profline8(4);  // T (4)
}


/* 16-bit version with 8 channels */

auto dprofile_shuffle16(WORD * dprofile,
                        WORD * score_matrix,
                        BYTE * dseq_byte) -> void
{
  static constexpr unsigned int channels {8};  // does 8 represent the number of channels?
  auto * profile_db = reinterpret_cast<__m128i *>(dprofile);
  auto * score_db = reinterpret_cast<__m128i *>(score_matrix);
  auto * sequence_db = reinterpret_cast<__m128i *>(dseq_byte);

  const auto zero = _mm_setzero_si128();
  const auto one = _mm_set1_epi16(1);

  auto transform_lower_seq_chunk = [&](const __m128i& seq_chunk) -> __m128i {
    auto lower_chunk = _mm_unpacklo_epi8(seq_chunk, zero);
    lower_chunk = _mm_slli_epi16(lower_chunk, 1);
    auto local_t = _mm_adds_epu16(lower_chunk, one);
    local_t = _mm_slli_epi16(local_t, channels);
    return _mm_or_si128(lower_chunk, local_t);
  };

  auto transform_higher_seq_chunk = [&](const __m128i& seq_chunk) -> __m128i {
    auto higher_chunk = _mm_unpackhi_epi8(seq_chunk, zero);
    higher_chunk = _mm_slli_epi16(higher_chunk, 1);
    auto local_t = _mm_adds_epu16(higher_chunk, one);
    local_t = _mm_slli_epi16(local_t, channels);
    return _mm_or_si128(higher_chunk, local_t);
  };

  const auto t0 = _mm_load_si128(sequence_db + 0);
  const auto m0 = transform_lower_seq_chunk(t0);
  const auto m1 = transform_higher_seq_chunk(t0);

  const auto t3 = _mm_load_si128(sequence_db + 1);
  const auto m2 = transform_lower_seq_chunk(t3);
  const auto m3 = transform_higher_seq_chunk(t3);

  auto profline16 = [&](const long long int nuc) {
    const auto scores = _mm_load_si128(score_db + 4 * nuc);

    _mm_store_si128(profile_db + 4 * nuc + 0, _mm_shuffle_epi8(scores, m0));
    _mm_store_si128(profile_db + 4 * nuc + 1, _mm_shuffle_epi8(scores, m1));
    _mm_store_si128(profile_db + 4 * nuc + 2, _mm_shuffle_epi8(scores, m2));
    _mm_store_si128(profile_db + 4 * nuc + 3, _mm_shuffle_epi8(scores, m3));
  };

  profline16(0);  // -/gap/no nucleotide (0)
  profline16(1);  // A (1)
  profline16(2);  // C (2)
  profline16(3);  // G (3)
  profline16(4);  // T (4)
}

#else
#error __SSSE3__ not defined
#endif
#endif
