/*
    SWARM

    Copyright (C) 2012-2016 Torbjorn Rognes and Frederic Mahe

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

#include <stdio.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <stdlib.h>
#include <regex.h>
#include <limits.h>
#include <city.h>
#include <sys/types.h>
#include <sys/resource.h>

#ifdef __APPLE__
#include <sys/sysctl.h>
#else
#include <sys/sysinfo.h>
#endif

#ifdef __SSE2__
#include <emmintrin.h>
#endif

#ifdef __SSSE3__
#include <tmmintrin.h>
#endif

/* constants */

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#define SWARM_VERSION "2.1.8"
#define WIDTH 32
#define WIDTH_SHIFT 5
#define BLOCKWIDTH 32
#define MAX_THREADS 256
#define SEPCHAR ' '

#ifdef BIASED
#define ZERO 0x00
#else
#define ZERO 0x80
#endif

#ifndef MIN
#define MIN(x,y) ((x)<(y)?(x):(y))
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
  unsigned int headerlen;
  unsigned int headeridlen;
  unsigned int seqlen;
  unsigned int abundance;
  unsigned int clusterid;
  unsigned int hdrhash;
  int abundance_start;
  int abundance_end;
};

typedef struct seqinfo_s seqinfo_t;

extern seqinfo_t * seqindex;
extern qgramvector_t * qgrams;

struct queryinfo
{
  unsigned long qno;
  long len;
  char * seq;
};

typedef struct queryinfo queryinfo_t;

/* common data */

extern char * queryname;
extern char * matrixname;
extern long gapopen;
extern long gapextend;
extern long gapopenextend;
extern long matchscore;
extern long mismatchscore;
extern unsigned long threads;
extern char * databasename;
extern long resolution;
extern long mothur;
extern long usearch_abundance;

extern char map_ncbi_nt4[];
extern char map_ncbi_nt16[];
extern char map_ncbi_aa[];
extern char map_sound[];

extern long penalty_factor;
extern long penalty_gapextend;
extern long penalty_gapopen;
extern long penalty_mismatch;

extern FILE * outfile;
extern FILE * statsfile;
extern FILE * uclustfile;
extern FILE * internal_structure_file;
extern FILE * logfile;
extern FILE * fp_seeds;

extern char * opt_log;
extern char * opt_internal_structure;
extern char * opt_seeds;
extern long opt_no_otu_breaking;
extern long opt_fastidious;
extern long opt_boundary;
extern long opt_bloom_bits;
extern long opt_ceiling;
extern long opt_append_abundance;

extern long SCORELIMIT_7;
extern long SCORELIMIT_8;
extern long SCORELIMIT_16;
extern long SCORELIMIT_32;
extern long SCORELIMIT_63;
extern char BIAS;

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;

extern unsigned char * score_matrix_8;
extern unsigned short * score_matrix_16;
extern long * score_matrix_63;

extern char sym_nt[];

extern unsigned long longestdbsequence;

extern queryinfo_t query;


/* functions in util.cc */

long gcd(long a, long b);
void fatal(const char * msg);
void fatal(const char * format, const char * message);
void * xmalloc(size_t size);
void * xrealloc(void * ptr, size_t size);
char * xstrchrnul(char *s, int c);
unsigned long hash_fnv_1a_64(unsigned char * s, unsigned long n);
unsigned int hash_fnv_1a_32(unsigned char * s, unsigned long n);
unsigned long hash_djb2(unsigned char * s, unsigned long n);
unsigned long hash_djb2a(unsigned char * s, unsigned long n);
unsigned long hash_cityhash64(unsigned char * s, unsigned long n);
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done();

/* functions in qgram.cc */

void findqgrams(unsigned char * seq, unsigned long seqlen,
                unsigned char * qgramvector);
unsigned long qgram_diff(unsigned long a, unsigned long b);
void qgram_diff_fast(unsigned long seed,
                     unsigned long listlen,
                     unsigned long * amplist,
                     unsigned long * difflist);
void qgram_diff_init();
void qgram_diff_done();

/* functions in db.cc */

void db_read(const char * filename);

unsigned long db_getsequencecount();
unsigned long db_getnucleotidecount();

unsigned long db_getlongestheader();
unsigned long db_getlongestsequence();

seqinfo_t * db_getseqinfo(unsigned long seqno);

char * db_getsequence(unsigned long seqno);
unsigned long db_getsequencelen(unsigned long seqno);

void db_getsequenceandlength(unsigned long seqno,
                             char ** address,
                             long * length);

char * db_getheader(unsigned long seqno);
unsigned long db_getheaderlen(unsigned long seqno);

unsigned long db_getabundance(unsigned long seqno);

void db_showsequence(unsigned long seqno);
void db_showall();
void db_free();

void db_putseq(long seqno);

void db_qgrams_init();
void db_qgrams_done();

void db_fprintseq(FILE * fp, int a, int width);

inline unsigned char * db_getqgramvector(unsigned long seqno)
{
  return (unsigned char*)(qgrams + seqno);
}

void fprint_id(FILE * stream, unsigned long x);
void fprint_id_noabundance(FILE * stream, unsigned long x);
void fprint_id_with_new_abundance(FILE * stream,
                                  unsigned long seqno,
                                  unsigned long abundance);


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
             unsigned long sequences,
             unsigned long * seqnos,
             unsigned long * scores,
             unsigned long * diffs,
             unsigned long * alignmentlengths,
             unsigned long qlen,
             unsigned long dirbuffersize,
             unsigned long * dirbuffer);


/* functions in search16.cc */

void search16(WORD * * q_start,
              WORD gap_open_penalty,
              WORD gap_extend_penalty,
              WORD * score_matrix,
              WORD * dprofile,
              WORD * hearray,
              unsigned long sequences,
              unsigned long * seqnos,
              unsigned long * scores,
              unsigned long * diffs,
              unsigned long * alignmentlengths,
              unsigned long qlen,
              unsigned long dirbuffersize,
              unsigned long * dirbuffer);


/* functions in nw.cc */

void nw(char * dseq,
        char * dend,
        char * qseq,
        char * qend,
        long * score_matrix,
        unsigned long gapopen,
        unsigned long gapextend,
        unsigned long * nwscore,
        unsigned long * nwdiff,
        unsigned long * nwalignmentlength,
        char ** nwalignment,
        unsigned char * dir,
        unsigned long * hearray,
        unsigned long queryno,
        unsigned long dbseqno);


/* functions in matrix.cc */

void score_matrix_init();
void score_matrix_free();


/* functions in scan.cc */

void search_all(unsigned long query_no);
void search_do(unsigned long query_no, 
               unsigned long listlength,
               unsigned long * targets,
               unsigned long * scores,
               unsigned long * diffs,
               unsigned long * alignlengths,
               long bits);
void search_begin();
void search_end();


/* functions in algo.cc */

void algo_run();
void algo_d1_run();


/* functions in derep.cc */

void dereplicate();


/* functions in arch.cc */

unsigned long arch_get_memused();
unsigned long arch_get_memtotal();


/* new header files */

#include "bitmap.h"
#include "bloom.h"
#include "threads.h"
