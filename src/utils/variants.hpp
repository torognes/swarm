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

#include "view.hpp"  // View
#include <cstdint>
#include <vector>


class Zobrist;       // defined in zobrist.hpp
struct Sequence;     // defined in db.hpp


/* Variant information */
enum struct Variant_type : unsigned char { substitution, deletion, insertion };

struct var_s
{
  uint64_t hash;
  unsigned int pos;
  Variant_type type;
  unsigned char base;  // encoded nucleotide: 0, 1, 2, or 3
  unsigned short dummy; /* for alignment padding only */
};

// The sequence of one microvariant of seed, written into the caller's
// buffer and returned as a view of it -- previously a length handed back
// through an out-parameter, which every caller had to pair with the
// buffer again to rebuild exactly this value. The buffer is sized once
// per thread (longest sequence plus one insertion) and reused for every
// variant, so the returned Sequence is valid only until the next call.
auto generate_variant_sequence(Sequence const & seed,
                               struct var_s const & var,
                               std::vector<char> & buffer) -> Sequence;

auto check_variant(Sequence const & seed,
                   struct var_s const & var,
                   Sequence const & amp) -> bool;

// The microvariants of seq, written into the caller's variant_list and
// returned as the window that was filled. variant_list is sized to the
// 7L+4 upper bound once per thread, so the return value is what says how
// much of it this call made valid -- previously a count that every caller
// had to pair with the buffer again.
auto generate_variants(Zobrist const & zobrist,
                       Sequence const & seq,
                       uint64_t hash,
                       std::vector<struct var_s> & variant_list) -> View<struct var_s>;
