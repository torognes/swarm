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

#ifndef SWARM_UTILS_QGRAM_ARRAY_H
#define SWARM_UTILS_QGRAM_ARRAY_H

#include <array>
#include <vector>


// 128 bytes = 1,024 bits, one bit per possible 5-mer (4^5 = 1,024).
// qgramvectorbytes is derived from qgramlength below.
constexpr unsigned int qgramlength     {5};
constexpr unsigned int qgramvectorbytes {(1U << (2 * qgramlength)) / 8};

using Qgram_vector = std::array<unsigned char, qgramvectorbytes>;
using Qgram_store  = std::vector<Qgram_vector>;

#endif
