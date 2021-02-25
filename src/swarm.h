/*
    SWARM

    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

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

#include <cinttypes>

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

#include <array>
#include <cassert>
#include <climits>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <getopt.h>
#include <pthread.h>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>  // replace with <fstream> to improve portability
#include <vector>

#ifdef __APPLE__
#include <sys/resource.h>
#include <sys/sysctl.h>
#elif defined _WIN32
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#include <sys/sysinfo.h>
#endif

#ifdef __aarch64__

#include <arm_neon.h>

#elif defined __x86_64__

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

#ifdef __SSE4_1__
#include <smmintrin.h>
#endif

#ifdef __POPCNT__
#include <popcntintrin.h>
#endif

#define CAST_m128i_ptr(x) (reinterpret_cast<__m128i*>(x))

#elif defined __PPC__

#ifdef __LITTLE_ENDIAN__
#include <altivec.h>
#else
#error Big endian ppc64 CPUs not supported
#endif

#else

#error Unknown architecture
#endif

static_assert(INT_MAX > INT16_MAX, "Your compiler uses very short integers.");

/* constants */

const std::string swarm_version = {"Swarm 3.0.0"};
// constexpr unsigned int width {32};  // unused?
// constexpr unsigned int width_shift {5};  // unused?
constexpr char sepchar {' '};
constexpr unsigned int bloom_bits_default {16};

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

constexpr unsigned int qgramlength {5};
constexpr unsigned int qgramvectorbits {1 << (2 * qgramlength)};
constexpr unsigned int qgramvectorbytes {qgramvectorbits / 8};

/* structures and data types */

using WORD = unsigned short;
using BYTE = unsigned char;

using qgramvector_t = unsigned char[qgramvectorbytes];

struct seqinfo_s
{
  char * header;
  char * seq;
  uint64_t abundance;
  uint64_t hdrhash;
  uint64_t seqhash;
  int headerlen;
  unsigned int seqlen;
  unsigned int clusterid;
  int abundance_start;
  int abundance_end;
  int dummy; /* alignment padding only */
};

using seqinfo_t = struct seqinfo_s;

extern seqinfo_t * seqindex;
extern qgramvector_t * qgrams;

struct queryinfo
{
  uint64_t qno;
  int64_t len;
  char * seq;
};

using queryinfo_t = struct queryinfo;

/* common data */

struct Parameters {
  bool opt_help {false};
  bool opt_disable_sse3 {false};
  bool opt_version {false};
  std::string input_filename {"-"};
  std::string opt_network_file;
  std::string opt_internal_structure;
  bool opt_fastidious {false};
  int64_t opt_bloom_bits {bloom_bits_default};
  int64_t opt_differences {1};
};

extern std::string opt_log;  // used by multithreaded functions
extern std::string opt_output_file;
extern std::string opt_seeds;
extern std::string opt_statistics_file;
extern std::string opt_uclust_file;
extern int64_t opt_append_abundance;
extern int64_t opt_boundary;  // used by multithreaded functions
extern int64_t opt_ceiling;
extern int64_t opt_gap_extension_penalty;
extern int64_t opt_gap_opening_penalty;
extern int64_t opt_match_reward;
extern int64_t opt_mismatch_penalty;
extern bool opt_mothur;
extern bool opt_no_otu_breaking;
extern int64_t opt_threads;
extern bool opt_usearch_abundance;

extern int64_t penalty_factor;
extern int64_t penalty_gapextend;
extern int64_t penalty_gapopen;
extern int64_t penalty_mismatch;

extern FILE * outfile;
extern FILE * statsfile;
extern FILE * uclustfile;
extern FILE * internal_structure_file;
extern FILE * logfile;
extern FILE * fp_seeds;
extern FILE * network_file;

extern int64_t mmx_present;
extern int64_t sse_present;
extern int64_t sse2_present;
extern int64_t sse3_present;
extern int64_t ssse3_present;
extern int64_t sse41_present;
extern int64_t sse42_present;
extern int64_t popcnt_present;
extern int64_t avx_present;
extern int64_t avx2_present;

extern unsigned char * score_matrix_8;
extern unsigned short * score_matrix_16;
extern int64_t * score_matrix_63;

extern const std::array<char, 32> sym_nt;

extern uint64_t longestdbsequence;

extern queryinfo_t query;

extern uint64_t duplicates_found;

/* inline functions */

inline auto nt_extract(char * seq, uint64_t i) -> unsigned char
{
  // Extract compressed nucleotide in sequence seq at position i
  constexpr unsigned int max_nt_per_uint64 {1 << 5};  // 32 nt fit in 64 bits
  constexpr unsigned int drop_remainder {5};  // (len+31) % 32 = remainder (drop it)
  constexpr unsigned int max_range {3};
  // outputs four possible values: 0, 1, 2 or 3
  return (((reinterpret_cast<uint64_t*>(seq))[i >> drop_remainder]) >> \
          ((i & (max_nt_per_uint64 - 1)) << 1)) & max_range;
}

inline auto nt_bytelength(unsigned int len) -> unsigned int
{
  // Compute number of bytes used for compressed sequence of length len
  // (minimum result is 8 bytes)
  constexpr unsigned int max_nt_per_uint64 {1 << 5};  // 32 nt fit in 64 bits
  constexpr unsigned int drop_remainder {5};  // (len+31) % 32 = remainder (drop it)
  constexpr unsigned int convert_to_bytes {3};  // times 8 to get the number of bytes
  return ((len + max_nt_per_uint64 - 1) >> drop_remainder) << convert_to_bytes;
}

/* functions in util.cc */

auto gcd(int64_t a, int64_t b) -> int64_t;
[[ noreturn ]] void fatal(const char * msg);
auto xmalloc(size_t size) -> void *;
auto xrealloc(void * ptr, size_t size) -> void *;
void xfree(void * ptr);
void progress_init(const char * prompt, uint64_t size);
void progress_update(uint64_t progress);
void progress_done();
auto fopen_input(const char * filename) -> FILE *;
auto fopen_output(const char * filename) -> FILE *;

/* functions in qgram.cc */

void findqgrams(unsigned char * seq, uint64_t seqlen,
                unsigned char * qgramvector);
auto qgram_diff(uint64_t a, uint64_t b) -> uint64_t;
void qgram_diff_fast(uint64_t seed,
                     uint64_t listlen,
                     uint64_t * amplist,
                     uint64_t * difflist);
void qgram_diff_init();
void qgram_diff_done();

/* functions in db.cc */

void db_read(const char * filename, struct Parameters const & p);

auto db_getsequencecount() -> unsigned int;
auto db_getnucleotidecount() -> uint64_t;

auto db_getlongestheader() -> unsigned int;
auto db_getlongestsequence() -> unsigned int;

auto db_getseqinfo(uint64_t seqno) -> seqinfo_t *;

auto db_getsequence(uint64_t seqno) -> char *;
auto db_getsequencelen(uint64_t seqno) -> unsigned int;

auto db_gethash(uint64_t seqno) -> uint64_t;

void db_getsequenceandlength(uint64_t seqno,
                             char ** address,
                             unsigned int * length);

auto db_getheader(uint64_t seqno) -> char *;
auto db_getheaderlen(uint64_t seqno) -> unsigned int;

auto db_getabundance(uint64_t seqno) -> uint64_t;

void db_showall();
void db_free();

void db_putseq(int64_t seqno);

void db_qgrams_init();
void db_qgrams_done();

void db_fprintseq(FILE * fp, unsigned int a, unsigned int width);

inline auto db_getqgramvector(uint64_t seqno) -> unsigned char *
{
  return reinterpret_cast<unsigned char*>(qgrams + seqno);
}

void fprint_id(FILE * stream, uint64_t x);
void fprint_id_noabundance(FILE * stream, uint64_t x);
void fprint_id_with_new_abundance(FILE * stream,
                                  uint64_t seqno,
                                  uint64_t abundance);

#ifdef __x86_64__

/* functions in ssse3.cc */

void dprofile_shuffle8(BYTE * dprofile,
                       BYTE * score_matrix,
                       BYTE * dseq_byte);

void dprofile_shuffle16(WORD * dprofile,
                        WORD * score_matrix,
                        BYTE * dseq_byte);

/* function in popcnt.cc */

auto compareqgramvectors_popcnt(unsigned char * a, unsigned char * b)
  -> uint64_t;

#endif

/* functions in search8.cc */

void search8(BYTE * * q_start,
             BYTE gap_open_penalty,
             BYTE gap_extend_penalty,
             BYTE * score_matrix,
             BYTE * dprofile,
             BYTE * hearray,
             uint64_t sequences,
             uint64_t * seqnos,
             uint64_t * scores,
             uint64_t * diffs,
             uint64_t * alignmentlengths,
             uint64_t qlen,
             uint64_t dirbuffersize,
             uint64_t * dirbuffer);


/* functions in search16.cc */

void search16(WORD * * q_start,
              WORD gap_open_penalty,
              WORD gap_extend_penalty,
              WORD * score_matrix,
              WORD * dprofile,
              WORD * hearray,
              uint64_t sequences,
              uint64_t * seqnos,
              uint64_t * scores,
              uint64_t * diffs,
              uint64_t * alignmentlengths,
              uint64_t qlen,
              uint64_t dirbuffersize,
              uint64_t * dirbuffer);


/* functions in nw.cc */

void nw(char * dseq,
        int64_t dlen,
        char * qseq,
        int64_t qlen,
        int64_t * score_matrix,
        int64_t gapopen,
        int64_t gapextend,
        int64_t * nwscore,
        int64_t * nwdiff,
        int64_t * nwalignmentlength,
        char ** nwalignment,
        unsigned char * dir,
        int64_t * hearray,
        uint64_t queryno,
        uint64_t dbseqno);


/* functions in matrix.cc */

void score_matrix_init();
void score_matrix_free();


/* functions in scan.cc */

void search_all(uint64_t query_no);
void search_do(uint64_t query_no,
               uint64_t listlength,
               uint64_t * targets,
               uint64_t * scores,
               uint64_t * diffs,
               uint64_t * alignlengths,
               int bits);
void search_begin();
void search_end();


/* functions in algo.cc */

void algo_run(struct Parameters const & p);
void algo_d1_run(struct Parameters const & p);


/* functions in derep.cc */

void dereplicate(struct Parameters const & p);


/* functions in arch.cc */

auto arch_get_memused() -> uint64_t;
auto arch_get_memtotal() -> uint64_t;
void arch_srandom(unsigned int seed);
auto arch_random() -> uint64_t;


/* new header files */

#include "threads.h"
#include "zobrist.h"
#include "bloompat.h"
#include "bloomflex.h"
#include "variants.h"
#include "hashtable.h"
