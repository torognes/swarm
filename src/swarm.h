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

#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // FILE
#include <string>


/* constants */

// C++17 refactor: use 'inline constexpr' in header files

// At file scope, constexpr implies const, and const implies
// static. So plain constexpr variables in headers will be duplicated
// in all cpp files that include them, and because they have internal
// linkage, they wont be de-duplicated.  'inline constexpr' variables
// will be de-duplicated by the linker.
// (source: https://www.youtube.com/watch?v=QVHwOOrSh3w)

constexpr char dash_filename {'-'};
constexpr unsigned int opt_differences_default {1};
constexpr unsigned int ceiling_default {0};
constexpr auto boundary_default = 3;
constexpr unsigned int append_abundance_default {0};
constexpr unsigned int mismatch_penalty_default {4};
constexpr unsigned int match_reward_default {5};
constexpr auto gap_opening_penalty_default = 12L;
constexpr unsigned int gap_extension_penalty_default {4};
constexpr unsigned int bloom_bits_default {16};


/* common data */

struct Parameters {
  int64_t opt_bloom_bits {bloom_bits_default};
  int64_t opt_differences {opt_differences_default};
  int64_t opt_mismatch_penalty {mismatch_penalty_default};
  int64_t opt_match_reward {match_reward_default};
  int64_t opt_gap_opening_penalty {gap_opening_penalty_default};
  int64_t opt_gap_extension_penalty {gap_extension_penalty_default};
  int64_t opt_ceiling {ceiling_default};
  int64_t opt_append_abundance {append_abundance_default};
  int64_t opt_boundary {boundary_default};
  int64_t mmx_present {0};
  int64_t sse42_present {0};
  int64_t sse_present {0};
  int64_t sse2_present {0};
  int64_t avx2_present {0};
  int64_t avx_present {0};
  int64_t sse3_present {0};
  int64_t sse41_present {0};
  int64_t ssse3_present {0};
  int64_t popcnt_present {0};
  int64_t penalty_mismatch {2 * match_reward_default +
                            2 * mismatch_penalty_default};
  int64_t penalty_gapextend {match_reward_default +
                             2 * gap_extension_penalty_default};
  int64_t penalty_gapopen {2 * gap_opening_penalty_default};
  bool opt_help {false};
  bool opt_disable_sse3 {false};
  bool opt_version {false};
  bool opt_fastidious {false};
  bool opt_usearch_abundance {false};
  bool opt_mothur {false};
  std::string input_filename {dash_filename};
  std::string opt_network_file;
  std::string opt_internal_structure;
  std::string opt_seeds;
  std::string opt_statistics_file;
  std::string opt_uclust_file;
  std::string opt_output_file {dash_filename};
};

// Note: extern - static storage duration and external linkage
//
// Objects of static storage duration are zero-initialized before
// main() is called (verified in gdb). However, it seems that the
// linker requires these objects to be explicitly initialized in one
// of the translation units.

extern std::string opt_log;  // used by multithreaded functions
extern bool opt_no_cluster_breaking;  // three function calls
extern int64_t opt_threads;

extern std::FILE * outfile;
extern std::FILE * statsfile;
extern std::FILE * uclustfile;
extern std::FILE * internal_structure_file;
extern std::FILE * logfile;
extern std::FILE * fp_seeds;
extern std::FILE * network_file;
