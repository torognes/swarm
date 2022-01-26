/*
    SWARM

    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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

#include <cstdint>
#include <smmintrin.h>


using VECTORTYPE = __m128i;


void align_cells_regular_16_sse41(VECTORTYPE * Sm,
                                  VECTORTYPE * hep,
                                  VECTORTYPE ** qp,
                                  VECTORTYPE * Qm,
                                  VECTORTYPE * Rm,
                                  uint64_t ql,
                                  VECTORTYPE * F0,
                                  uint64_t * dir_long,
                                  VECTORTYPE * H0);

void align_cells_masked_16_sse41(VECTORTYPE * Sm,
                                 VECTORTYPE * hep,
                                 VECTORTYPE ** qp,
                                 VECTORTYPE * Qm,
                                 VECTORTYPE * Rm,
                                 uint64_t ql,
                                 VECTORTYPE * F0,
                                 uint64_t * dir_long,
                                 VECTORTYPE * H0,
                                 VECTORTYPE * Mm,
                                 VECTORTYPE * MQ,
                                 VECTORTYPE * MR,
                                 VECTORTYPE * MQ0);
