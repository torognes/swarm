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

#include "utils/seqinfo.h"
#include "utils/view.h"
#include "utils/zobrist.h"
#include <cstdio>  // std::FILE
#include <cstdint>  // uint64_t
#include <memory>  // std::unique_ptr
#include <vector>


struct Parameters;  // defined in swarm.h


// Non-owning view of a packed-nucleotide amplicon. `length` is the
// nucleotide count (canonical: matches the historical seqlen).
// `encoded` views the storage bytes, with encoded.size() ==
// nt_bytelength(length) (4 nt per byte). Downstream consumers
// reinterpret encoded.data() as uint64_t* for SIMD-friendly access;
// iterating encoded directly walks the packed bytes, not nucleotides.
struct Sequence {
  View<char> encoded;
  unsigned int length;
};


class Data {
public:
  explicit Data(struct Parameters const & parameters);

  auto sequence_count()   const -> unsigned int { return static_cast<unsigned int>(seqindex_.size()); }
  auto longest_sequence() const -> unsigned int { return longest_; }

  auto zobrist() const -> Zobrist const & { return *zobrist_p_; }

  auto info(uint64_t seqno) const -> struct seqinfo_s const &;

  auto sequence_view(uint64_t seqno)   const -> Sequence;
  auto sequence_hash(uint64_t seqno)   const -> uint64_t;
  auto header_view(uint64_t seqno)     const -> View<char>;
  auto abundance(uint64_t seqno)       const -> uint64_t;

  auto fprintseq(std::FILE * stream, unsigned int seqno) const -> void;
  auto fprint_id(std::FILE * stream,
                 uint64_t seqno,
                 bool opt_usearch_abundance,
                 int64_t opt_append_abundance) const -> void;
  auto fprint_id_noabundance(std::FILE * stream,
                             uint64_t seqno,
                             bool opt_usearch_abundance) const -> void;
  auto fprint_id_with_new_abundance(std::FILE * stream,
                                    uint64_t seqno,
                                    uint64_t new_abundance,
                                    bool opt_usearch_abundance) const -> void;

private:
  std::vector<char>             data_;
  std::vector<struct seqinfo_s> seqindex_;
  std::unique_ptr<Zobrist>      zobrist_p_;  // deferred: needs longest_sequence
  unsigned int                  longest_ {0};
};


