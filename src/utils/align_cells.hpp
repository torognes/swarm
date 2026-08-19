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


#ifndef SWARM_UTILS_ALIGN_CELLS_H
#define SWARM_UTILS_ALIGN_CELLS_H


// The alignment kernel, in one place: one body for all three of the builds
// that used to hold a copy of it.
//
// They differed in two compile-time facts and nothing else. The width -- an
// 8-bit lane with 16 channels (search8.cpp) against a 16-bit lane with 8
// (search16.cpp) -- and, at 16 bits on x86-64, whether an unsigned minimum is
// one PMINUW or has to be emulated (arch/x86_64/sse41.cpp against the
// baseline). Both now arrive as one policy type, 'Ops', the way the
// masked/regular flavour already arrives as the type of 'masks'.
//
// Ops supplies, as static members:
//
//   Dir_word            the direction-word type: uint16_t at both widths,
//                       because a mask carries one bit per channel and there
//                       are never more than 16
//   add, sub, min       lane arithmetic
//   mask_eq             lanewise compare, packed into a Dir_word
//   zero                an all-zero vector
//
// A class rather than function pointers, so the calls disappear at
// instantiation -- checked in the disassembly, not assumed.
//
// This header names, but deliberately does not include:
//
//   VECTORTYPE     the channel type, whose spelling is per-architecture *and*
//                  per-width (__m128i, uint8x16_t, uint16x8_t, __vector
//                  unsigned char, __vector unsigned short) and which must not
//                  become a template parameter -- naming it as one makes GCC
//                  drop its attributes, may_alias among them (see
//                  utils/mask_vectors.hpp);
//   No_mask        from utils/mask_vectors.hpp: the one piece of the masking
//                  machinery that is genuinely width-independent;
//   max_ptrdiff    used by the assertions, and only defined in debug builds.
//
// So it has to be included after those are in scope, which is why it carries
// no includes of its own and why each includer names it last, inside its own
// anonymous namespace -- which also gives the instantiations the internal
// linkage the hand-written copies had.


// The masking payload: 'mask' selects the channels whose sequence just
// ended, 'mq' is the running gap-open accumulator seeded by the caller,
// 'mr' its per-iteration increment, and 'mq0' the value 'mq' held on
// entry. A plain struct rather than a template on VECTORTYPE, which cannot be
// a template argument here at all -- see utils/mask_vectors.hpp.
struct Mask_vectors {
  VECTORTYPE mask;
  VECTORTYPE mq;
  VECTORTYPE mr;
  VECTORTYPE mq0;
};


// The masking step, selected by the type of the kernel's mask argument
// (see utils/mask_vectors.hpp). The No_mask overload is empty, so the
// regular kernel's loop body contains nothing at this point. Both are
// templates on the width policy because the masked one arithmetics with it;
// the call site names it explicitly, and the mask parameter still chooses
// between them.
template <typename Ops>
inline auto apply_mask(VECTORTYPE & /*h4*/, VECTORTYPE & /*E*/,
                       No_mask const & /*masks*/) -> void
{
}

template <typename Ops>
inline auto apply_mask(VECTORTYPE & h4, VECTORTYPE & E,
                       Mask_vectors & masks) -> void
{
  /* mask h4 and E */
  h4 = Ops::sub(h4, masks.mask);
  E  = Ops::sub(E,  masks.mask);

  /* init h4 and E */
  h4 = Ops::add(h4, masks.mq);
  E  = Ops::add(E,  masks.mq);
  E  = Ops::add(E,  masks.mq0);

  /* update MQ */
  masks.mq = Ops::add(masks.mq,  masks.mr);
}


// One cell of the block. The DIR writes use a subscript rather than
// std::next for the reason given on the loop below.
template <typename Ops>
inline auto onestep(VECTORTYPE & H,
                       VECTORTYPE & N,
                       VECTORTYPE & F,
                       VECTORTYPE const V,
                       typename Ops::Dir_word * const DIR,
                       VECTORTYPE & E,
                       VECTORTYPE const QR,
                       VECTORTYPE const R) -> void
{
  H = Ops::add(H, V);
  auto const W = H;
  H = Ops::min(H, F);
  DIR[0] = Ops::mask_eq(W, H);  // subscript, not std::next: hot loop, see align_cells
  H = Ops::min(H, E);
  DIR[1] = Ops::mask_eq(H, E);
  N = H;
  H = Ops::add(H, QR);
  F = Ops::add(F, R);
  E = Ops::add(E, R);
  F = Ops::min(H, F);
  DIR[2] = Ops::mask_eq(H, F);
  E = Ops::min(H, E);
  DIR[3] = Ops::mask_eq(H, E);
}


// One block of cells, shared by the regular and masked kernels. The
// masked variant differs only by a per-iteration adjustment of h4 and E;
// which flavour this is comes from the type of 'masks', so the regular
// instantiation drops that adjustment entirely and is handed no masking
// data at all (see utils/mask_vectors.hpp).
template <typename Ops, typename Masks>
auto align_cells(VECTORTYPE * const Sm,
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

  auto * const dir = reinterpret_cast<typename Ops::Dir_word *>(dir_long);

  auto const Q = Qm;
  auto const R = Rm;

  auto f0 = F0;
  auto f1 = Ops::add(f0, R);
  auto f2 = Ops::add(f1, R);
  auto f3 = Ops::add(f2, R);

  auto h0 = H0;
  auto h1 = Ops::sub(f0, Q);
  auto h2 = Ops::add(h1, R);
  auto h3 = Ops::add(h2, R);

  auto h5 = Ops::zero();
  auto h6 = Ops::zero();
  auto h7 = Ops::zero();
  auto h8 = Ops::zero();

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

      apply_mask<Ops>(h4, E, masks);

      onestep<Ops>(h0, h5, f0, x[0], &dir[(step * pos) + offset0], E, Q, R);
      onestep<Ops>(h1, h6, f1, x[1], &dir[(step * pos) + offset1], E, Q, R);
      onestep<Ops>(h2, h7, f2, x[2], &dir[(step * pos) + offset2], E, Q, R);
      onestep<Ops>(h3, h8, f3, x[3], &dir[(step * pos) + offset3], E, Q, R);
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

#endif  // SWARM_UTILS_ALIGN_CELLS_H
