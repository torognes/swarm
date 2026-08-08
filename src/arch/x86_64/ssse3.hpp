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

#ifndef SWARM_ARCH_X86_64_SSSE3_H
#define SWARM_ARCH_X86_64_SSSE3_H


#ifdef __SSSE3__

#include "../../utils/search_data.hpp"  // BYTE, WORD, Score_matrix_8/16, Dseq_8/16

auto dprofile_shuffle8(Dprofile_8 & dprofile_a,
                       Score_matrix_8 const & score_matrix_a,
                       Dseq_8 const & dseq_a) -> void;

auto dprofile_shuffle16(Dprofile_16 & dprofile_a,
                        Score_matrix_16 const & score_matrix_a,
                        Dseq_16 const & dseq_a) -> void;

#endif

#endif  // SWARM_ARCH_X86_64_SSSE3_H
