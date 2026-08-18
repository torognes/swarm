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


#ifndef SWARM_UTILS_ALIGN_CELLS_16_H
#define SWARM_UTILS_ALIGN_CELLS_16_H


// The 16-bit alignment kernel, in one place.
//
// It exists in two builds that differ by a single instruction: the
// baseline uses v_min16, which on x86-64 without SSE4.1 has to emulate an
// unsigned 16-bit minimum, and the SSE4.1 build uses _mm_min_epu16
// (PMINUW). Which one applies is a compile-time fact -- it is the
// instruction set the translation unit is built for -- so it arrives as a
// template argument, the way the masked/regular flavour already arrives as
// the type of 'masks' (see mask_vectors.hpp). Before this, the fact was
// expressed by keeping a second copy of the whole ~490-token body in
// arch/x86_64/sse41.cpp, and the two copies had begun to drift.
//
// Min_op is a class with one static member, min(lhs, rhs); a class rather
// than a function pointer or a std::function so that the call disappears
// at instantiation.
//
// This header names, but deliberately does not include:
//
//   VECTORTYPE     the channel type, whose spelling is per-architecture
//                  (__m128i, uint16x8_t, __vector unsigned short) and
//                  which must not become a template parameter -- naming it
//                  as one makes GCC drop its attributes, may_alias among
//                  them (see utils/mask_vectors.hpp);
//   WORD           the direction-word type;
//   v_add16, v_sub16, v_zero16, v_mask_eq16
//                  the intrinsic wrappers, which are distinct names per
//                  width rather than overloads;
//   Mask_vectors, apply_mask
//                  per-kernel by design, again see mask_vectors.hpp;
//   max_ptrdiff    used by the assertions, and only defined in debug builds.
//
// So it has to be included after all of those are in scope, which is why
// it carries no includes of its own and why the two includers name it
// last. That is the price of not templating on the vector type.


// One cell of the block. The DIR writes use a subscript rather than
// std::next for the reason given on the loop below.
template <typename Min_op>
inline auto onestep_16(VECTORTYPE & H,
                       VECTORTYPE & N,
                       VECTORTYPE & F,
                       VECTORTYPE const V,
                       WORD * const DIR,
                       VECTORTYPE & E,
                       VECTORTYPE const QR,
                       VECTORTYPE const R) -> void
{
  H = v_add16(H, V);
  auto const W = H;
  H = Min_op::min(H, F);
  DIR[0] = v_mask_eq16(W, H);  // subscript, not std::next: hot loop, see align_cells
  H = Min_op::min(H, E);
  DIR[1] = v_mask_eq16(H, E);
  N = H;
  H = v_add16(H, QR);
  F = v_add16(F, R);
  E = v_add16(E, R);
  F = Min_op::min(H, F);
  DIR[2] = v_mask_eq16(H, F);
  E = Min_op::min(H, E);
  DIR[3] = v_mask_eq16(H, E);
}


// One block of cells, shared by the regular and masked kernels. The
// masked variant differs only by a per-iteration adjustment of h4 and E;
// which flavour this is comes from the type of 'masks', so the regular
// instantiation drops that adjustment entirely and is handed no masking
// data at all (see utils/mask_vectors.hpp).
template <typename Min_op, typename Masks>
auto align_cells_16(VECTORTYPE * const Sm,
                    VECTORTYPE * const hep,
                    VECTORTYPE ** const qp,
                    VECTORTYPE const & Qm,
                    VECTORTYPE const & Rm,
                    uint64_t const ql,
                    VECTORTYPE const & F0,
                    uint64_t * const dir_long,
                    VECTORTYPE const & H0,
                    Masks & masks) -> void
{
  static constexpr auto step = 16;
  static constexpr auto offset0 = 0;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  VECTORTYPE E;
  VECTORTYPE h4;

  auto * const dir = reinterpret_cast<WORD *>(dir_long);

  auto const Q = Qm;
  auto const R = Rm;

  auto f0 = F0;
  auto f1 = v_add16(f0, R);
  auto f2 = v_add16(f1, R);
  auto f3 = v_add16(f2, R);

  auto h0 = H0;
  auto h1 = v_sub16(f0, Q);
  auto h2 = v_add16(h1, R);
  auto h3 = v_add16(h2, R);

  auto h5 = v_zero16();
  auto h6 = v_zero16();
  auto h7 = v_zero16();
  auto h8 = v_zero16();

  assert(ql <= max_ptrdiff);
  assert(ql <= ((max_ptrdiff - 1) / 2));  // max 'E' offset
  assert(ql <= ((max_ptrdiff - offset3) / step));  // max 'dir' offset
  auto const ql_signed = static_cast<std::ptrdiff_t>(ql);
  // Performance: subscript / &dir[...] rather than std::next() in this hot
  // loop. The std::next() form (commit 8c6925f, taken for
  // cppcoreguidelines-pro-bounds-pointer-arithmetic) pessimized the SSE4.1
  // kernel by ~20 % on d > 1 18SV9. Stays clang-tidy clean anyway: pos is
  // signed, so no -Wsign-conversion, and operator[] is not pointer
  // arithmetic.
  for (auto pos = 0LL; pos < ql_signed; ++pos)
    {
      VECTORTYPE const * const x = qp[pos];
      h4 = hep[(2 * pos) + 0];
      E  = hep[(2 * pos) + 1];

      apply_mask(h4, E, masks);

      onestep_16<Min_op>(h0, h5, f0, x[0], &dir[(step * pos) + offset0], E, Q, R);
      onestep_16<Min_op>(h1, h6, f1, x[1], &dir[(step * pos) + offset1], E, Q, R);
      onestep_16<Min_op>(h2, h7, f2, x[2], &dir[(step * pos) + offset2], E, Q, R);
      onestep_16<Min_op>(h3, h8, f3, x[3], &dir[(step * pos) + offset3], E, Q, R);
      hep[(2 * pos) + 0] = h8;
      hep[(2 * pos) + 1] = E;
      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}

#endif  // SWARM_UTILS_ALIGN_CELLS_16_H
