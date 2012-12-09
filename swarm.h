/*
    SWARM

    Copyright (C) 2012 Torbjorn Rognes and Frederic Mahe

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
#include <tmmintrin.h>
#include <stdlib.h>

/* constants */

#ifndef LINE_MAX
#define LINE_MAX 2048
#endif

#define SWARM_VERSION "1.0.2"
#define WIDTH 32
#define WIDTH_SHIFT 5
#define BLOCKWIDTH 32
#define MAX_THREADS 256

#ifdef BIASED
#define ZERO 0x00
#else
#define ZERO 0x80
#endif

#define MIN(x,y) ((x)<(y)?(x):(y))

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;
typedef BYTE VECTOR[16];

struct seqinfo_s
{
  char * header;
  char * seq;
  unsigned long headerlen;
  unsigned long seqlen;
  unsigned long abundance;
  unsigned long clusterid;
};

typedef struct seqinfo_s seqinfo_t;

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
extern long threads;
extern char * databasename;
extern long resolution;

extern char map_ncbi_nt4[];
extern char map_ncbi_nt16[];
extern char map_ncbi_aa[];
extern char map_sound[];

extern long penalty_factor;
extern long penalty_gapextend;
extern long penalty_gapopen;
extern long penalty_mismatch;

extern FILE * out;

extern long SCORELIMIT_7;
extern long SCORELIMIT_8;
extern long SCORELIMIT_16;
extern long SCORELIMIT_32;
extern long SCORELIMIT_63;
extern char BIAS;

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

/* functions in dataset.cc */

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

/* functions in search8.cc, search16.cc, nw.cc */

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
             unsigned long qlen,
	     unsigned long dirbuffersize,
	     WORD * up,
	     WORD * left);

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
	      unsigned long qlen,
	      unsigned long dirbuffersize,
	      WORD * up,
	      WORD * left);

void nw(char * dseq,
	char * dend,
	char * qseq,
	char * qend,
	unsigned long * hearray,
	unsigned long * score_matrix,
	unsigned long gap_open_penalty,
	unsigned long gap_extend_penalty,
	unsigned long * nwscore,
	unsigned long * nwdiff);

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
	       long bits);
void search_begin();
void search_end();

/* functions in algo.cc */

void algo_run();
