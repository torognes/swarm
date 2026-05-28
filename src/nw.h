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

#include <array>
#include <cstdint>  // int64_t
#include <string>
#include <vector>


// refactor: 'n_cells' is already defined in 'score_matrix.h'
constexpr auto n_cells_ = 32ULL;  // number of chars in sym_nt

auto nw(char const * dseq,
        uint64_t dlen,
        char const * qseq,
        uint64_t qlen,
        const std::array<int64_t, n_cells_ * n_cells_> & score_matrix,
        uint64_t gapopen,
        uint64_t gapextend,
        uint64_t & nwdiff,
        std::vector<unsigned char> & directions,
        std::vector<uint64_t> & hearray,
        std::vector<char> & raw_alignment) -> void;


/*
  Needleman/Wunsch/Sellers aligner

  Owns the scratch buffers (directions, hearray, raw_alignment) and the
  cigar output buffer, so callers don't have to allocate or clear them
  between calls. Sized once at construction with the longest sequence
  length; subsequent align() calls reuse the storage.

  align() returns a small Result aggregating the four values consumers
  need: cigar string, number of differences, alignment length, and
  percent identity. The cigar_string reference lives in this Alignment
  object's internal buffer and is only valid until the next align()
  call on the same object.
*/
class Alignment {
public:
  struct Result {
    std::string const & cigar_string;
    uint64_t            differences;
    uint64_t            length;       // == nwalignmentlength
    double              percent_id;
  };

  explicit Alignment(uint64_t longest_sequence);

  auto align(char const * dseq, uint64_t dlen,
             char const * qseq, uint64_t qlen,
             std::array<int64_t, n_cells_ * n_cells_> const & score_matrix,
             uint64_t gapopen,
             uint64_t gapextend) -> Result;

private:
  std::vector<unsigned char> directions_;
  std::vector<uint64_t>      hearray_;
  std::vector<char>          raw_alignment_;
  std::string                cigar_string_;
};
