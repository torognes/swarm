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

#include <array>
#include <cstdint>  // int64_t


// refactor: 'n_cells' is already defined in 'score_matrix.h'
constexpr auto n_cells_ = 32ULL;  // number of chars in sym_nt

void nw(char * dseq,
        int64_t dlen,
        char * qseq,
        int64_t qlen,
        const std::array<int64_t, n_cells_ * n_cells_> & score_matrix,
        int64_t gapopen,
        int64_t gapextend,
        int64_t * nwscore,
        int64_t * nwdiff,
        int64_t * nwalignmentlength,
        char ** nwalignment,
        unsigned char * dir,
        int64_t * hearray,
        uint64_t queryno,
        uint64_t dbseqno);
