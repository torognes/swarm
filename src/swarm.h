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

#include <cinttypes>  // macros PRIu64 and PRId64

#ifndef PRIu64
#ifdef _WIN32
#define PRIu64 "I64u"
#else
constexpr char PRIu64[] = "lu";
#endif
#endif

#ifndef PRId64
#ifdef _WIN32
#define PRId64 "I64d"
#else
constexpr char PRId64[] = "ld";
#endif
#endif

#include <cstdarg>  // refactoring: C-style variadic functions (unused?)
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <unistd.h>  // replace with <fstream> to improve portability

#ifdef __aarch64__
#include <arm_neon.h>
#elif defined __x86_64__

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#else

#error Unknown architecture
#endif

/* constants */

// C++17 refactor: use 'inline constexpr' in header files

// At file scope, constexpr implies const, and const implies
// static. So plain constexpr variables in headers will be duplicated
// in all cpp files that include them, and because they have internal
// linkage, they wont be de-duplicated.  'inline constexpr' variables
// will be de-duplicated by the linker.
// (source: https://www.youtube.com/watch?v=QVHwOOrSh3w)

const std::string swarm_version {"3.1.3"};
constexpr char sepchar {' '};
constexpr char dash_filename {'-'};
constexpr unsigned int opt_differences_default {1};
constexpr unsigned int ceiling_default {0};
constexpr unsigned int append_abundance_default {0};
constexpr unsigned int mismatch_penalty_default {4};
constexpr unsigned int match_reward_default {5};
constexpr unsigned int gap_opening_penalty_default {12};
constexpr unsigned int gap_extension_penalty_default {4};
constexpr unsigned int bloom_bits_default {16};
constexpr unsigned int qgramlength {5};
constexpr unsigned int qgramvectorbits {1U << (2 * qgramlength)};
constexpr unsigned int qgramvectorbytes {qgramvectorbits / 8};


/* structures and data types */

using qgramvector_t = unsigned char[qgramvectorbytes];
extern qgramvector_t * qgrams;

struct queryinfo
{
  uint64_t qno;
  int64_t len;
  char * seq;
};

using queryinfo_t = struct queryinfo;
extern queryinfo_t query;

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
  int64_t mmx_present {0};
  int64_t sse42_present {0};
  int64_t sse_present {0};
  int64_t sse2_present {0};
  int64_t avx2_present {0};
  int64_t avx_present {0};
  int64_t sse3_present {0};
  int64_t penalty_mismatch {2 * match_reward_default + 2 * mismatch_penalty_default};
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

extern std::string opt_log;  // used by multithreaded functions
extern int64_t opt_boundary;  // used by multithreaded functions
extern bool opt_no_otu_breaking;  // three function calls
extern int64_t opt_threads;

extern int64_t penalty_gapextend;
extern int64_t penalty_gapopen;
// extern int64_t penalty_mismatch;

extern std::FILE * outfile;
extern std::FILE * statsfile;
extern std::FILE * uclustfile;
extern std::FILE * internal_structure_file;
extern std::FILE * logfile;
extern std::FILE * fp_seeds;
extern std::FILE * network_file;

extern int64_t ssse3_present; // several function calls
extern int64_t sse41_present; // several function calls
extern int64_t popcnt_present; // several function calls

extern unsigned char * score_matrix_8;
extern unsigned short * score_matrix_16;
extern int64_t * score_matrix_63;
