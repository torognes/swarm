/*
    SWARM

    Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe

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

#include <inttypes.h>

#ifndef PRIu64
#ifdef _WIN32
#define PRIu64 "I64u"
#else
#define PRIu64 "lu"
#endif
#endif

#ifndef PRId64
#ifdef _WIN32
#define PRId64 "I64d"
#else
#define PRId64 "ld"
#endif
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <stdlib.h>
#include <regex.h>
#include <limits.h>
#include <stdarg.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "city.h"

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

static_assert(INT_MAX > 32767, "Your compiler uses very short integers.");

/* constants */

#define SWARM_VERSION "3.0.0"
#define WIDTH 32
#define WIDTH_SHIFT 5
#define BLOCKWIDTH 32
const unsigned int MAX_THREADS = 256;
#define SEPCHAR ' '

#ifdef BIASED
#define ZERO 0x00
#else
#define ZERO 0x80
#endif

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
#endif

#ifndef MAX
#define MAX(x,y) ((x)>(y)?(x):(y))
#endif

#define QGRAMLENGTH 5
#define QGRAMVECTORBITS (1<<(2*QGRAMLENGTH))
#define QGRAMVECTORBYTES (QGRAMVECTORBITS/8)

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

typedef unsigned char qgramvector_t[QGRAMVECTORBYTES];

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

typedef struct seqinfo_s seqinfo_t;

extern seqinfo_t * seqindex;
extern qgramvector_t * qgrams;

struct queryinfo
{
  uint64_t qno;
  int64_t len;
  char * seq;
};

typedef struct queryinfo queryinfo_t;

/* common data */

extern char * opt_internal_structure;
extern char * opt_log;
extern char * opt_output_file;
extern char * opt_seeds;
extern char * opt_statistics_file;
extern char * opt_uclust_file;
extern int64_t opt_append_abundance;
extern int64_t opt_bloom_bits;
extern int64_t opt_boundary;
extern int64_t opt_ceiling;
extern int64_t opt_differences;
extern int64_t opt_fastidious;
extern int64_t opt_gap_extension_penalty;
extern int64_t opt_gap_opening_penalty;
extern int64_t opt_help;
extern int64_t opt_match_reward;
extern int64_t opt_mismatch_penalty;
extern int64_t opt_mothur;
extern int64_t opt_no_otu_breaking;
extern int64_t opt_threads;
extern int64_t opt_usearch_abundance;
extern int64_t opt_version;

extern char * queryname;
extern char * matrixname;

extern char * input_filename;

extern char map_ncbi_nt4[];
extern char map_ncbi_nt16[];
extern char map_ncbi_aa[];
extern char map_sound[];

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


extern int64_t SCORELIMIT_7;
extern int64_t SCORELIMIT_8;
extern int64_t SCORELIMIT_16;
extern int64_t SCORELIMIT_32;
extern int64_t SCORELIMIT_63;
extern char BIAS;

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

extern char sym_nt[];

extern uint64_t longestdbsequence;

extern queryinfo_t query;

extern uint64_t duplicates_found;

/* inline functions */

inline unsigned char nt_extract(char * seq, uint64_t i)
{
  // Extract compressed nucleotide in sequence seq at position i
  return (((reinterpret_cast<uint64_t*>(seq))[i >> 5]) >> ((i & 31) << 1)) & 3;
}

inline unsigned int nt_bytelength(unsigned int len)
{
  // Compute number of bytes used for compressed sequence of length len
  return ((len+31) >> 5) << 3;
}

/* functions in util.cc */

int64_t gcd(int64_t a, int64_t b);
[[ noreturn ]] void fatal(const char * msg);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
void xfree(void * ptr);
uint64_t hash_fnv_1a_64(unsigned char * s, uint64_t n);
unsigned int hash_fnv_1a_32(unsigned char * s, uint64_t n);
uint64_t hash_djb2(unsigned char * s, uint64_t n);
uint64_t hash_djb2a(unsigned char * s, uint64_t n);
uint64_t hash_cityhash64(unsigned char * s, uint64_t n);
uint64_t hash_xor64len(unsigned char * s, uint64_t n);
uint64_t hash64shift(uint64_t key);
void progress_init(const char * prompt, uint64_t size);
void progress_update(uint64_t progress);
void progress_done();
FILE * fopen_input(const char * filename);
FILE * fopen_output(const char * filename);

/* functions in qgram.cc */

void findqgrams(unsigned char * seq, uint64_t seqlen,
                unsigned char * qgramvector);
uint64_t qgram_diff(uint64_t a, uint64_t b);
void qgram_diff_fast(uint64_t seed,
                     uint64_t listlen,
                     uint64_t * amplist,
                     uint64_t * difflist);
void qgram_diff_init();
void qgram_diff_done();

/* functions in db.cc */

void db_read(const char * filename);

unsigned int db_getsequencecount();
uint64_t db_getnucleotidecount();

unsigned int db_getlongestheader();
unsigned int db_getlongestsequence();

seqinfo_t * db_getseqinfo(uint64_t seqno);

char * db_getsequence(uint64_t seqno);
unsigned int db_getsequencelen(uint64_t seqno);

uint64_t db_gethash(uint64_t seqno);

void db_getsequenceandlength(uint64_t seqno,
                             char ** address,
                             unsigned int * length);

char * db_getheader(uint64_t seqno);
unsigned int db_getheaderlen(uint64_t seqno);

uint64_t db_getabundance(uint64_t seqno);

void db_showall();
void db_free();

void db_putseq(int64_t seqno);

void db_qgrams_init();
void db_qgrams_done();

void db_fprintseq(FILE * fp, unsigned int a, unsigned int width);

inline unsigned char * db_getqgramvector(uint64_t seqno)
{
  return reinterpret_cast<unsigned char*>(qgrams + seqno);
}

void fprint_id(FILE * stream, uint64_t x);
void fprint_id_noabundance(FILE * stream, uint64_t x);
void fprint_id_with_new_abundance(FILE * stream,
                                  uint64_t seqno,
                                  uint64_t abundance);


/* functions in ssse3.cc */

void dprofile_shuffle8(BYTE * dprofile,
                       BYTE * score_matrix,
                       BYTE * dseq_byte);

void dprofile_shuffle16(WORD * dprofile,
                        WORD * score_matrix,
                        BYTE * dseq_byte);


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

void algo_run();
void algo_d1_run();


/* functions in derep.cc */

void dereplicate();


/* functions in arch.cc */

uint64_t arch_get_memused();
uint64_t arch_get_memtotal();
void arch_srandom(unsigned int seed);
uint64_t arch_random();


/* new header files */

#include "threads.h"
#include "zobrist.h"
#include "bloompat.h"
#include "bloomflex.h"
#include "variants.h"
#include "hashtable.h"
