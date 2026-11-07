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

#include "align_rvv.hpp"

#ifdef __riscv_vector

#include "../../utils/simd_alignment.hpp"  // simd_vector_bytes
#include "intrinsics_to_functions.hpp"  // v_u16_t, v_u8_t
#include <riscv_vector.h>
#include <cassert>
#include <cstddef>  // std::ptrdiff_t, std::size_t
#include <cstdint>  // uint16_t, uint64_t, uint8_t
#ifndef NDEBUG
#include <limits>
#endif


// The kernel of utils/align_cells.hpp, written with RVV intrinsics: one
// instruction each for the saturating add/sub, and each direction mask
// stored straight from the mask register (vsm), where the generic version
// needs a reduction per mask. The RVV types are sizeless, so the vectors
// live in registers here and the buffers keep their generic types.
//
// vl is fixed at 128 bits, the width the buffers are laid out for (16 or 8
// channels); the V extension guarantees VLEN >= 128.
//
// Built with -fno-lto (see the Makefile): GCC 14 crashes on the mask
// intrinsics under LTO.

namespace {

#ifndef NDEBUG
constexpr auto max_ptrdiff = std::numeric_limits<std::ptrdiff_t>::max();
#endif

// direction word, as in utils/align_cells.hpp
using Dir_word = unsigned short;


template <typename Lane>
struct Rvv;

template <>
struct Rvv<uint8_t> {
  using Generic = v_u8_t;
  using Vec = vuint8m1_t;
  static constexpr std::size_t lanes {simd_vector_bytes / sizeof(uint8_t)};
  static auto load(Generic const * const src) -> Vec {
    return __riscv_vle8_v_u8m1(reinterpret_cast<uint8_t const *>(src), lanes);
  }
  static auto store(Generic * const dst, Vec const vec) -> void {
    __riscv_vse8_v_u8m1(reinterpret_cast<uint8_t *>(dst), vec, lanes);
  }
  static auto zero() -> Vec { return __riscv_vmv_v_x_u8m1(0, lanes); }
};

template <>
struct Rvv<uint16_t> {
  using Generic = v_u16_t;
  using Vec = vuint16m1_t;
  static constexpr std::size_t lanes {simd_vector_bytes / sizeof(uint16_t)};
  static auto load(Generic const * const src) -> Vec {
    return __riscv_vle16_v_u16m1(reinterpret_cast<uint16_t const *>(src), lanes);
  }
  static auto store(Generic * const dst, Vec const vec) -> void {
    __riscv_vse16_v_u16m1(reinterpret_cast<uint16_t *>(dst), vec, lanes);
  }
  static auto zero() -> Vec { return __riscv_vmv_v_x_u16m1(0, lanes); }
};

template <typename Lane>
using Vec = typename Rvv<Lane>::Vec;

template <typename Lane>
using Generic = typename Rvv<Lane>::Generic;


// One bit per channel, channel 0 in bit 0: 16 bits at 8-bit lanes, the low
// 8 bits of the word at 16-bit lanes (the high byte is not written).
template <typename Mask>
inline auto store_mask(Dir_word * const dir, Mask const mask,
                       std::size_t const n_lanes) -> void
{
  __riscv_vsm(reinterpret_cast<uint8_t *>(dir), mask, n_lanes);
}


template <typename Lane>
inline auto onestep(Vec<Lane> & H,
                    Vec<Lane> & N,
                    Vec<Lane> & F,
                    Vec<Lane> const V,
                    Dir_word * const DIR,
                    Vec<Lane> & E,
                    Vec<Lane> const QR,
                    Vec<Lane> const R) -> void
{
  static constexpr auto n_lanes = Rvv<Lane>::lanes;
  H = __riscv_vsaddu(H, V, n_lanes);
  auto const W = H;
  H = __riscv_vminu(H, F, n_lanes);
  store_mask(&DIR[0], __riscv_vmseq(W, H, n_lanes), n_lanes);
  H = __riscv_vminu(H, E, n_lanes);
  store_mask(&DIR[1], __riscv_vmseq(H, E, n_lanes), n_lanes);
  N = H;
  H = __riscv_vsaddu(H, QR, n_lanes);
  F = __riscv_vsaddu(F, R, n_lanes);
  E = __riscv_vsaddu(E, R, n_lanes);
  F = __riscv_vminu(H, F, n_lanes);
  store_mask(&DIR[2], __riscv_vmseq(H, F, n_lanes), n_lanes);
  E = __riscv_vminu(H, E, n_lanes);
  store_mask(&DIR[3], __riscv_vmseq(H, E, n_lanes), n_lanes);
}


// The regular and masked kernels; 'masked' is a compile-time constant, so
// the regular instantiation drops the masking step and its unused vectors.
template <typename Lane, bool masked>
inline auto align_rvv(Generic<Lane> * const Sm,
                      Generic<Lane> * const hep,
                      Generic<Lane> ** const qp,
                      Generic<Lane> const & Qm,
                      Generic<Lane> const & Rm,
                      uint64_t const ql,
                      Generic<Lane> const & F0,
                      uint64_t * const dir_long,
                      Generic<Lane> const & H0,
                      Vec<Lane> const mask,
                      Vec<Lane> & mq,
                      Vec<Lane> const mr,
                      Vec<Lane> const mq0) -> void
{
  using Ops = Rvv<Lane>;
  static constexpr auto n_lanes = Ops::lanes;
  static constexpr auto step = 16;
  static constexpr auto offset0 = 0;
  static constexpr auto offset1 = offset0 + 4;
  static constexpr auto offset2 = offset1 + 4;
  static constexpr auto offset3 = offset2 + 4;

  auto * const dir = reinterpret_cast<Dir_word *>(dir_long);

  auto const Q = Ops::load(&Qm);
  auto const R = Ops::load(&Rm);

  auto f0 = Ops::load(&F0);
  auto f1 = __riscv_vsaddu(f0, R, n_lanes);
  auto f2 = __riscv_vsaddu(f1, R, n_lanes);
  auto f3 = __riscv_vsaddu(f2, R, n_lanes);

  auto h0 = Ops::load(&H0);
  auto h1 = __riscv_vssubu(f0, Q, n_lanes);
  auto h2 = __riscv_vsaddu(h1, R, n_lanes);
  auto h3 = __riscv_vsaddu(h2, R, n_lanes);

  auto h5 = Ops::zero();
  auto h6 = Ops::zero();
  auto h7 = Ops::zero();
  auto h8 = Ops::zero();

  assert(ql <= max_ptrdiff);
  assert(ql <= ((max_ptrdiff - 1) / 2));  // max 'E' offset
  assert(ql <= ((max_ptrdiff - offset3) / step));  // max 'dir' offset
  auto const ql_signed = static_cast<std::ptrdiff_t>(ql);
  for (auto pos = 0LL; pos < ql_signed; ++pos)
    {
      Generic<Lane> const * const x = qp[pos];
      auto h4 = Ops::load(&hep[(2 * pos) + 0]);
      auto E  = Ops::load(&hep[(2 * pos) + 1]);

      if (masked) {
        h4 = __riscv_vssubu(h4, mask, n_lanes);
        E  = __riscv_vssubu(E,  mask, n_lanes);
        h4 = __riscv_vsaddu(h4, mq, n_lanes);
        E  = __riscv_vsaddu(E,  mq, n_lanes);
        E  = __riscv_vsaddu(E,  mq0, n_lanes);
        mq = __riscv_vsaddu(mq, mr, n_lanes);
      }

      onestep<Lane>(h0, h5, f0, Ops::load(&x[0]), &dir[(step * pos) + offset0], E, Q, R);
      onestep<Lane>(h1, h6, f1, Ops::load(&x[1]), &dir[(step * pos) + offset1], E, Q, R);
      onestep<Lane>(h2, h7, f2, Ops::load(&x[2]), &dir[(step * pos) + offset2], E, Q, R);
      onestep<Lane>(h3, h8, f3, Ops::load(&x[3]), &dir[(step * pos) + offset3], E, Q, R);
      Ops::store(&hep[(2 * pos) + 0], h8);
      Ops::store(&hep[(2 * pos) + 1], E);
      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  Ops::store(&Sm[0], h5);
  Ops::store(&Sm[1], h6);
  Ops::store(&Sm[2], h7);
  Ops::store(&Sm[3], h8);
}

}  // namespace


auto align_rvv_regular_8(v_u8_t * const Sm, v_u8_t * const hep,
                         v_u8_t ** const qp,
                         v_u8_t const & Qm, v_u8_t const & Rm,
                         uint64_t const ql, v_u8_t const & F0,
                         uint64_t * const dir_long, v_u8_t const & H0) -> void
{
  using Ops = Rvv<uint8_t>;
  auto mq = Ops::zero();
  align_rvv<uint8_t, false>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0,
                            Ops::zero(), mq, Ops::zero(), Ops::zero());
}


auto align_rvv_masked_8(v_u8_t * const Sm, v_u8_t * const hep,
                        v_u8_t ** const qp,
                        v_u8_t const & Qm, v_u8_t const & Rm,
                        uint64_t const ql, v_u8_t const & F0,
                        uint64_t * const dir_long, v_u8_t const & H0,
                        v_u8_t const * const Mm, v_u8_t * const MQ,
                        v_u8_t const * const MR, v_u8_t const * const MQ0) -> void
{
  using Ops = Rvv<uint8_t>;
  auto mq = Ops::load(MQ);
  align_rvv<uint8_t, true>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0,
                           Ops::load(Mm), mq, Ops::load(MR), Ops::load(MQ0));
  Ops::store(MQ, mq);
}


auto align_rvv_regular_16(v_u16_t * const Sm, v_u16_t * const hep,
                          v_u16_t ** const qp,
                          v_u16_t const & Qm, v_u16_t const & Rm,
                          uint64_t const ql, v_u16_t const & F0,
                          uint64_t * const dir_long, v_u16_t const & H0) -> void
{
  using Ops = Rvv<uint16_t>;
  auto mq = Ops::zero();
  align_rvv<uint16_t, false>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0,
                             Ops::zero(), mq, Ops::zero(), Ops::zero());
}


auto align_rvv_masked_16(v_u16_t * const Sm, v_u16_t * const hep,
                         v_u16_t ** const qp,
                         v_u16_t const & Qm, v_u16_t const & Rm,
                         uint64_t const ql, v_u16_t const & F0,
                         uint64_t * const dir_long, v_u16_t const & H0,
                         v_u16_t const * const Mm, v_u16_t * const MQ,
                         v_u16_t const * const MR, v_u16_t const * const MQ0) -> void
{
  using Ops = Rvv<uint16_t>;
  auto mq = Ops::load(MQ);
  align_rvv<uint16_t, true>(Sm, hep, qp, Qm, Rm, ql, F0, dir_long, H0,
                            Ops::load(Mm), mq, Ops::load(MR), Ops::load(MQ0));
  Ops::store(MQ, mq);
}

#endif  // __riscv_vector
