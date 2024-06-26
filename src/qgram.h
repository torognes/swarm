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

#include <cstdint>  // uint64_t
#include <vector>


auto findqgrams(unsigned char * seq, uint64_t seqlen,
                unsigned char * qgramvector) -> void;
auto qgram_diff(uint64_t seqno_a, uint64_t seqno_b) -> uint64_t;
auto qgram_diff_fast(uint64_t seed,
                     uint64_t listlen,
                     uint64_t * amplist,
                     uint64_t * difflist,
                     std::vector<struct thread_info_s>& thread_info_v) -> void;
auto qgram_diff_init(std::vector<struct thread_info_s>& thread_info_v) -> void;
auto qgram_diff_done() -> void;
