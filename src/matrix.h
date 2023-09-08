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

#include <algorithm>
#include <array>
#include <cstdint>  // int64_t, uint16_t
#include <type_traits>

constexpr auto n_cells = 32ULL;  // number of chars in sym_nt

// note: there is no uchar8_t, only char8_t in C++20
// refactoring: C++20 use 'requires' to constrain accepted types
template <typename Integral>
auto create_score_matrix(const std::int64_t penalty_mismatch)
  -> std::array<Integral, n_cells * n_cells> {
  static_assert(std::is_same<Integral, unsigned char>::value \
                || std::is_same<Integral, unsigned short int>::value \
                || std::is_same<Integral, std::uint16_t>::value \
                || std::is_same<Integral, std::int64_t>::value,
                "Invalid type! Only unsigned char, unsigned short and int64_t can be used.");
  static constexpr Integral matchscore {0};
  const auto mismatchscore = static_cast<Integral>(penalty_mismatch);
  std::array<Integral, n_cells * n_cells> score_matrix {{}};

  // fill in the upper-left quarter of the matrix with mismatchscore,
  // except the diagonal starting from cell (1, 1)
  auto index = 0U;
  auto choose_score = [&index, &mismatchscore](Integral &element) {
    const auto column = index % n_cells;
    const auto row = index / n_cells;
    element = ((row == column && row != 0 && column != 0)
               || (column >= n_cells / 2)
               || (row >= n_cells / 2)) ? matchscore : mismatchscore;
    ++index;
  };
  // refactoring: C++20 ranges::for_each
  std::for_each(score_matrix.begin(), score_matrix.end(), choose_score);

  return score_matrix;
}


// usage example:
//
// create_score_matrix<unsigned char>(penalty_mismatch)};   OK
// create_score_matrix<unsigned short>(penalty_mismatch)};  OK
// create_score_matrix<int64_t>(penalty_mismatch)};         OK
// create_score_matrix<signed char>(penalty_mismatch)};     compilation error

// refactoring: C++20 constexpr template
// static_assert(create_score_matrix<unsigned short>(4)[0] == 4);  // (0, 0)
// static_assert(create_score_matrix<unsigned short>(4)[(n_cells / 2) - 1] == 4);  // (0, 15)
// static_assert(create_score_matrix<unsigned short>(4)[(n_cells / 2)] == 0);  // (0, 16)
// static_assert(create_score_matrix<unsigned short>(4)[(1 * n_cells) + 0] == 4);  // (1, 0)
// static_assert(create_score_matrix<unsigned short>(4)[(1 * n_cells) + 1] == 0);  // (1, 1)
// static_assert(create_score_matrix<unsigned short>(4)[(1 * n_cells) + 2] == 4);  // (1, 2)
// static_assert(create_score_matrix<unsigned short>(4)[(1 * n_cells) + (n_cells / 2) - 1] == 4);  // (1, 15)
// static_assert(create_score_matrix<unsigned short>(4)[(1 * n_cells) + (n_cells / 2)] == 0);  // (1, 16)
// static_assert(create_score_matrix<unsigned short>(4)[(((n_cells / 2) - 1) * n_cells) + (n_cells / 2) - 1] == 0);  // (15, 15)
// static_assert(create_score_matrix<unsigned short>(4)[(n_cells * n_cells) - 1] == 0);  // last cell


// expected score matrix (if mismatch score is 4, and type is unsigned short):
//
//    0 . . . . 5 . . . .10 . . . .15 . . . .20 . . . .25 . . . .30 .
// 0  4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 1  4 0 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 2  4 4 0 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 3  4 4 4 0 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 4  4 4 4 4 0 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 5  4 4 4 4 4 0 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 6  4 4 4 4 4 4 0 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 7  4 4 4 4 4 4 4 0 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 8  4 4 4 4 4 4 4 4 0 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 9  4 4 4 4 4 4 4 4 4 0 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 10 4 4 4 4 4 4 4 4 4 4 0 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 11 4 4 4 4 4 4 4 4 4 4 4 0 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 12 4 4 4 4 4 4 4 4 4 4 4 4 0 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 13 4 4 4 4 4 4 4 4 4 4 4 4 4 0 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 14 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 15 4 4 4 4 4 4 4 4 4 4 4 4 4 4 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 16 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 17 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 18 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 19 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 21 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 22 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 23 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 24 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 25 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 26 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 27 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 28 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 29 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 30 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
// 31 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
