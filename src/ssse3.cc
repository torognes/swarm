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


#ifdef __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>  // SSE2 intrinsics
#include "utils/intrinsics_to_functions_x86_64.h"
#include <iterator>  // std::next
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

auto v_shuffle8(__m128i lhs, __m128i mask) -> __m128i {
  // shuffle vector of bytes according to control mask
  return _mm_shuffle_epi8(lhs, mask);
}

auto v_or(__m128i lhs, __m128i rhs) -> __m128i {
  return _mm_or_si128(lhs, rhs);
}

auto v_shift_left(__m128i vector, int n_bytes) -> __m128i {
  // shift vector of shorts to the left by n bytes, pad with zeros
  return _mm_slli_epi16(vector, n_bytes);
}

// refactoring: why not using that in vzero8 and vzero16?
auto v_zero() -> __m128i {
  // return vector with all elements set to zero
  return _mm_setzero_si128();
}


/* 8-bit version with 16 channels */

auto dprofile_shuffle8(BYTE * dprofile,
                       BYTE * score_matrix,
                       BYTE * dseq_byte) -> void
{
  // inputs: score_matrix and dseq_byte (sequence from db); output: dprofile
  auto * const sequence_db = cast_vector8(dseq_byte);
  auto * const score_db = cast_vector8(score_matrix);
  auto * const profile_db = cast_vector8(dprofile);    // output
  const auto seq_chunk0 = v_load8(std::next(sequence_db, 0));  // 16 nucleotides
  const auto seq_chunk1 = v_load8(std::next(sequence_db, 1));  // next 16
  const auto seq_chunk2 = v_load8(std::next(sequence_db, 2));  // next 16
  const auto seq_chunk3 = v_load8(std::next(sequence_db, 3));  // final 16 (total of 64)

  auto profline8 = [&](const long long int nuc) {
    // scores: 16 scores from the score matrix, matching the
    // nucleotide 'nuc'; five different nucleotides (0, 1, 2, 3, 4),
    // so five possible rows of scores
    const auto scores = v_load8(std::next(score_db, 2 * nuc));

    v_store8(std::next(profile_db, (4 * nuc) + 0), v_shuffle8(scores, seq_chunk0));
    v_store8(std::next(profile_db, (4 * nuc) + 1), v_shuffle8(scores, seq_chunk1));
    v_store8(std::next(profile_db, (4 * nuc) + 2), v_shuffle8(scores, seq_chunk2));
    v_store8(std::next(profile_db, (4 * nuc) + 3), v_shuffle8(scores, seq_chunk3));
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
  // inputs: score_matrix and dseq_byte (sequence from db); output: dprofile
  auto * const profile_db = cast_vector16(dprofile);
  auto * const score_db = cast_vector16(score_matrix);
  auto * const sequence_db = cast_vector8(dseq_byte);
  static constexpr int channels {8};  // does 8 represent the number of channels?

  const auto zero = v_zero();
  const auto one = v_dup16(1);

  // refactoring: make lower_chunk and local_t const?
  auto transform_lower_seq_chunk = [&](const __m128i& seq_chunk) -> __m128i {
    auto lower_chunk = v_merge_lo_8(seq_chunk, zero);
    lower_chunk = v_shift_left(lower_chunk, 1);
    auto local_t = v_add16(lower_chunk, one);
    local_t = v_shift_left(local_t, channels);
    return v_or(lower_chunk, local_t);
  };

  auto transform_higher_seq_chunk = [&](const __m128i& seq_chunk) -> __m128i {
    auto higher_chunk = v_merge_hi_8(seq_chunk, zero);
    higher_chunk = v_shift_left(higher_chunk, 1);
    auto local_t = v_add16(higher_chunk, one);
    local_t = v_shift_left(local_t, channels);
    return v_or(higher_chunk, local_t);
  };

  const auto t0 = v_load8(std::next(sequence_db, 0));
  const auto m0 = transform_lower_seq_chunk(t0);
  const auto m1 = transform_higher_seq_chunk(t0);

  const auto t3 = v_load8(std::next(sequence_db, 1));
  const auto m2 = transform_lower_seq_chunk(t3);
  const auto m3 = transform_higher_seq_chunk(t3);

  auto profline16 = [&](const long long int nuc) {
    const auto scores = v_load16(std::next(score_db, 4 * nuc));

    v_store16(std::next(profile_db, (4 * nuc) + 0), v_shuffle8(scores, m0));
    v_store16(std::next(profile_db, (4 * nuc) + 1), v_shuffle8(scores, m1));
    v_store16(std::next(profile_db, (4 * nuc) + 2), v_shuffle8(scores, m2));
    v_store16(std::next(profile_db, (4 * nuc) + 3), v_shuffle8(scores, m3));
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
