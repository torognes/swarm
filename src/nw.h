/*
    SWARM

    Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

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
#include <vector>


// refactor: 'n_cells' is already defined in 'score_matrix.h'
constexpr auto n_cells_ = 32ULL;  // number of chars in sym_nt

auto nw(char * dseq,
        uint64_t dlen,
        char * qseq,
        uint64_t qlen,
        const std::array<int64_t, n_cells_ * n_cells_> & score_matrix,
        uint64_t gapopen,
        uint64_t gapextend,
        uint64_t & nwdiff,
        std::vector<unsigned char> & directions,
        std::vector<uint64_t> & hearray,
        std::vector<char> & raw_alignment) -> void;
