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

#ifndef SWARM_UTILS_MASK_VECTORS_H
#define SWARM_UTILS_MASK_VECTORS_H


// Mask state for the alignment kernels.
//
// The kernels come in two flavours: a regular one, and a masked one that
// adjusts h4 and E once per iteration when a database sequence ended in
// the middle of a block. Both share a single loop body, so the loop has
// to be told which flavour it is running.
//
// The flavour is carried by the *type* of the mask argument rather than
// by a value: the regular kernel is instantiated with No_mask, whose
// apply_mask() overload is empty, and the masking step vanishes. That
// keeps the two kernels sharing one body without a compile-time bool to
// thread through, and without handing the regular kernel four null
// pointers (or four dummy zero vectors) that it must not look at.
//
// Only the empty tag lives here, because it is the only part that is
// genuinely common. Two things stay next to each kernel:
//
//   - apply_mask(), which is written in terms of the width-specific
//     intrinsic wrappers (v_add8 / v_sub8 against v_add16 / v_sub16).
//     Those are distinct names rather than overloads, so they cannot be
//     selected generically;
//
//   - Mask_vectors, the four-vector payload, which cannot be a template
//     on the channel type. Naming VECTORTYPE as a template argument
//     makes GCC drop its attributes ("ignoring attributes on template
//     argument 'VECTORTYPE'"), and one of those attributes is
//     may_alias. This codebase already hit that wall: see the C array
//     'VECTORTYPE S[4]' in search8.cpp / search16.cpp, kept as a C array
//     for exactly this reason. Since these vectors sit next to the
//     H0 / F0 lane writes whose aliasing behaviour is delicate (see
//     set_lane_16), a plain per-kernel struct is the safe form.

struct No_mask { };

#endif
