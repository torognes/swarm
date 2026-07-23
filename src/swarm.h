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

#include "utils/input_output.h"  // FileHandle
#include <cstdint>  // int64_t, uint64_t
#include <cstdio>  // FILE, stderr
#include <string>


/* common data */

struct Parameters {
  // Defaults are scoped to Parameters: every reference is from inside
  // the struct, so an inner home avoids polluting the global namespace
  // and side-steps the file-scope constexpr ODR question entirely
  // (static constexpr class members have no duplication issue).
  static constexpr char dash_filename {'-'};
  static constexpr unsigned int opt_differences_default {1};
  static constexpr unsigned int ceiling_default {0};
  static constexpr auto boundary_default = 3;
  static constexpr unsigned int append_abundance_default {0};
  static constexpr unsigned int mismatch_penalty_default {4};
  static constexpr unsigned int match_reward_default {5};
  static constexpr auto gap_opening_penalty_default = 12L;
  static constexpr unsigned int gap_extension_penalty_default {4};
  static constexpr unsigned int bloom_bits_default {16};

  std::uint32_t opt_threads {1};
  std::uint64_t opt_bloom_bits {bloom_bits_default};
  std::uint64_t opt_differences {opt_differences_default};
  int64_t opt_mismatch_penalty {mismatch_penalty_default};
  int64_t opt_match_reward {match_reward_default};
  int64_t opt_gap_opening_penalty {gap_opening_penalty_default};
  int64_t opt_gap_extension_penalty {gap_extension_penalty_default};
  std::uint64_t opt_ceiling {ceiling_default};
  int64_t opt_append_abundance {append_abundance_default};
  std::uint64_t opt_boundary {boundary_default};
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
  int64_t penalty_mismatch {(2 * match_reward_default) + (2 * mismatch_penalty_default)};
  int64_t penalty_gapextend {match_reward_default + (2 * gap_extension_penalty_default)};
  int64_t penalty_gapopen {2 * gap_opening_penalty_default};
  bool opt_help {false};
  bool opt_disable_sse3 {false};
  bool opt_version {false};
  bool opt_fastidious {false};
  bool opt_usearch_abundance {false};
  bool opt_mothur {false};
  bool opt_no_cluster_breaking {false};
  std::string input_filename {dash_filename};
  std::string opt_network_file;
  std::string opt_internal_structure;
  std::string opt_seeds;
  std::string opt_statistics_file;
  std::string opt_uclust_file;
  std::string opt_output_file {dash_filename};
  std::string opt_log;
  // Output files. The FileHandle destructor closes the underlying
  // std::FILE * at scope exit; call sites pass <name>.get() to fprintf
  // and friends. logfile is the exception: it defaults to stderr and
  // is only optionally backed by an owned FileHandle (logfile_handle),
  // so it stays a raw pointer.
  FileHandle outfile;
  FileHandle statsfile;
  FileHandle uclustfile;
  FileHandle internal_structure_file;
  FileHandle seeds_file;
  FileHandle network_file;
  FileHandle logfile_handle;  // empty unless -l was given
  std::FILE * logfile {stderr};  // stderr macro expands to type std::FILE*
};
