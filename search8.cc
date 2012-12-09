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

#include "swarm.h"

// #define DEBUG

#define CHANNELS 16
#define CDEPTH 4

#define MATRIXWIDTH 16

void dprofile_dump8(BYTE * dprofile)
{
  char * ss = sym_nt;
  //  char * ss = sym_sound;

  printf("\ndprofile:\n");
  for(int k=0; k<4; k++)
  {
    printf("k=%d 0 1 2 3 4 5 6 7 8 9 a b c d e f\n", k);
    for(int i=0; i<16; i++)
    {
      printf("%c: ",ss[i]);
      for(int j=0; j<16; j++)
	printf("%2d", (char) dprofile[i*64+16*k+j]);
      printf("\n");
    }
  }
  printf("\n");
  exit(1);
}

int dumpcounter = 0;
char lines[4*16*1000];

void dseq_dump8(BYTE * dseq)
{
  char * s = sym_nt;

  if (dumpcounter < 21)
  {
    for(int i=0; i<CHANNELS; i++)
    {
      for(int j=0; j<CDEPTH; j++)
      {
	lines[4000*i+4*dumpcounter+j] = s[dseq[j*CHANNELS+i]];
      }
    }
    dumpcounter++;
  }
  else
  {
    for(int i=0; i<16; i++)
    {
      printf("%.1000s\n", lines+4000*i);
    }
    exit(1);
  }
}


inline void dprofile_shuffle8(BYTE * dprofile,
			      BYTE * score_matrix,
			      BYTE * dseq_byte)
{
  __m128i m0, m1, m2, m3, t0, t1, t2, t3, t4;

  __m128i * dseq = (__m128i*) dseq_byte;
  
  // ca 12 * 5 + 4 = 64 instructions

  // 16 x 4 = 64 db symbols

  // make masks

  /* Note: pshufb only on modern Intel cpus (SSSE3), not AMD */
  /* SSSE3: Supplemental SSE3 */

  m0 = _mm_load_si128(dseq);
  m1 = _mm_load_si128(dseq+1);
  m2 = _mm_load_si128(dseq+2);
  m3 = _mm_load_si128(dseq+3);

#define profline(j)					\
  t0 = _mm_load_si128((__m128i*)(score_matrix)+2*j);	\
  t1 = _mm_shuffle_epi8(t0, m0);			\
  t2 = _mm_shuffle_epi8(t0, m1);			\
  t3 = _mm_shuffle_epi8(t0, m2);			\
  t4 = _mm_shuffle_epi8(t0, m3);			\
  _mm_store_si128((__m128i*)(dprofile)+4*j+0, t1);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+1, t2);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+2, t3);	\
  _mm_store_si128((__m128i*)(dprofile)+4*j+3, t4)


  profline(0);
  profline(1);
  profline(2);
  profline(3);
  profline(4);
  /*
  profline(5);
  profline(6);
  profline(7);
  profline(8);
  profline(9);
  profline(10);
  profline(11);
  profline(12);
  profline(13);
  profline(14);
  profline(15);
  */

  //  dprofile_dump7(dprofile);
}

inline void dprofile_fill8(BYTE * dprofile,
			   BYTE * score_matrix,
			   BYTE * dseq)
{
  __m128i xmm0,  xmm1, xmm2,  xmm3,  xmm4,  xmm5,  xmm6,  xmm7;
  __m128i xmm8,  xmm9, xmm10, xmm11, xmm12, xmm13, xmm14, xmm15;
  
  // 4 x 16 db symbols
  // ca (60x2+68x2)x4 = 976 instructions

  for(int j=0; j<CDEPTH; j++)
  {
    unsigned d[CHANNELS];
    for(int i=0; i<CHANNELS; i++)
      d[i] = dseq[j*CHANNELS+i] << 5;
      
    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + d[0] ));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + d[2] ));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + d[4] ));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + d[6] ));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + d[8] ));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + d[10]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + d[12]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + d[14]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  *(__m128i*)(score_matrix + d[1] ));
    xmm2  = _mm_unpacklo_epi8(xmm2,  *(__m128i*)(score_matrix + d[3] ));
    xmm4  = _mm_unpacklo_epi8(xmm4,  *(__m128i*)(score_matrix + d[5] ));
    xmm6  = _mm_unpacklo_epi8(xmm6,  *(__m128i*)(score_matrix + d[7] ));
    xmm8  = _mm_unpacklo_epi8(xmm8,  *(__m128i*)(score_matrix + d[9] ));
    xmm10 = _mm_unpacklo_epi8(xmm10, *(__m128i*)(score_matrix + d[11]));
    xmm12 = _mm_unpacklo_epi8(xmm12, *(__m128i*)(score_matrix + d[13]));
    xmm14 = _mm_unpacklo_epi8(xmm14, *(__m128i*)(score_matrix + d[15]));
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128((__m128i*)(dprofile+16*j+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+448), xmm15);


    // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[0 ]));
    xmm1  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[1 ]));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[2 ]));
    xmm3  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[3 ]));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[4 ]));
    xmm5  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[5 ]));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[6 ]));
    xmm7  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[7 ]));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[8 ]));
    xmm9  = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[9 ]));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[10]));
    xmm11 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[11]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[12]));
    xmm13 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[13]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[14]));
    xmm15 = _mm_loadl_epi64((__m128i*)(score_matrix + 8 + d[15]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  xmm1);
    xmm2  = _mm_unpacklo_epi8(xmm2,  xmm3);
    xmm4  = _mm_unpacklo_epi8(xmm4,  xmm5);
    xmm6  = _mm_unpacklo_epi8(xmm6,  xmm7);
    xmm8  = _mm_unpacklo_epi8(xmm8,  xmm9);
    xmm10 = _mm_unpacklo_epi8(xmm10, xmm11);
    xmm12 = _mm_unpacklo_epi8(xmm12, xmm13);
    xmm14 = _mm_unpacklo_epi8(xmm14, xmm15);
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128((__m128i*)(dprofile+16*j+512+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+512+448), xmm15);


    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[0 ]));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[2 ]));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[4 ]));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[6 ]));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[8 ]));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[10]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[12]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + 16 + d[14]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  *(__m128i*)(score_matrix + 16 + d[1 ]));
    xmm2  = _mm_unpacklo_epi8(xmm2,  *(__m128i*)(score_matrix + 16 + d[3 ]));
    xmm4  = _mm_unpacklo_epi8(xmm4,  *(__m128i*)(score_matrix + 16 + d[5 ]));
    xmm6  = _mm_unpacklo_epi8(xmm6,  *(__m128i*)(score_matrix + 16 + d[7 ]));
    xmm8  = _mm_unpacklo_epi8(xmm8,  *(__m128i*)(score_matrix + 16 + d[9 ]));
    xmm10 = _mm_unpacklo_epi8(xmm10, *(__m128i*)(score_matrix + 16 + d[11 ]));
    xmm12 = _mm_unpacklo_epi8(xmm12, *(__m128i*)(score_matrix + 16 + d[13 ]));
    xmm14 = _mm_unpacklo_epi8(xmm14, *(__m128i*)(score_matrix + 16 + d[15 ]));
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128((__m128i*)(dprofile+16*j+1024+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+1024+448), xmm15);


    // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

    xmm0  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[0 ]));
    xmm1  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[1 ]));
    xmm2  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[2 ]));
    xmm3  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[3 ]));
    xmm4  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[4 ]));
    xmm5  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[5 ]));
    xmm6  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[6 ]));
    xmm7  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[7 ]));
    xmm8  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[8 ]));
    xmm9  = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[9 ]));
    xmm10 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[10]));
    xmm11 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[11]));
    xmm12 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[12]));
    xmm13 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[13]));
    xmm14 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[14]));
    xmm15 = _mm_loadl_epi64((__m128i*)(score_matrix + 24 + d[15]));

    xmm0  = _mm_unpacklo_epi8(xmm0,  xmm1);
    xmm2  = _mm_unpacklo_epi8(xmm2,  xmm3);
    xmm4  = _mm_unpacklo_epi8(xmm4,  xmm5);
    xmm6  = _mm_unpacklo_epi8(xmm6,  xmm7);
    xmm8  = _mm_unpacklo_epi8(xmm8,  xmm9);
    xmm10 = _mm_unpacklo_epi8(xmm10, xmm11);
    xmm12 = _mm_unpacklo_epi8(xmm12, xmm13);
    xmm14 = _mm_unpacklo_epi8(xmm14, xmm15);
      
    xmm1 = xmm0;
    xmm0 = _mm_unpacklo_epi16(xmm0, xmm2);
    xmm1 = _mm_unpackhi_epi16(xmm1, xmm2);
    xmm5 = xmm4;
    xmm4 = _mm_unpacklo_epi16(xmm4, xmm6);
    xmm5 = _mm_unpackhi_epi16(xmm5, xmm6);
    xmm9 = xmm8;
    xmm8 = _mm_unpacklo_epi16(xmm8, xmm10);
    xmm9 = _mm_unpackhi_epi16(xmm9, xmm10);
    xmm13 = xmm12;
    xmm12 = _mm_unpacklo_epi16(xmm12, xmm14);
    xmm13 = _mm_unpackhi_epi16(xmm13, xmm14);

    xmm2  = xmm0;
    xmm0  = _mm_unpacklo_epi32(xmm0, xmm4);
    xmm2  = _mm_unpackhi_epi32(xmm2, xmm4);
    xmm6  = xmm1;
    xmm1  = _mm_unpacklo_epi32(xmm1, xmm5);
    xmm6  = _mm_unpackhi_epi32(xmm6, xmm5);
    xmm10 = xmm8;
    xmm8  = _mm_unpacklo_epi32(xmm8, xmm12);
    xmm10 = _mm_unpackhi_epi32(xmm10, xmm12);
    xmm14 = xmm9;
    xmm9  = _mm_unpacklo_epi32(xmm9, xmm13);
    xmm14 = _mm_unpackhi_epi32(xmm14, xmm13);
      
    xmm3  = xmm0;
    xmm0  = _mm_unpacklo_epi64(xmm0, xmm8);
    xmm3  = _mm_unpackhi_epi64(xmm3, xmm8);
    xmm7  = xmm2;
    xmm2  = _mm_unpacklo_epi64(xmm2, xmm10);
    xmm7  = _mm_unpackhi_epi64(xmm7, xmm10);
    xmm11 = xmm1;
    xmm1  = _mm_unpacklo_epi64(xmm1, xmm9);
    xmm11 = _mm_unpackhi_epi64(xmm11, xmm9);
    xmm15 = xmm6;
    xmm6  = _mm_unpacklo_epi64(xmm6, xmm14);
    xmm15 = _mm_unpackhi_epi64(xmm15, xmm14);

    _mm_store_si128((__m128i*)(dprofile+16*j+1536+  0), xmm0);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+ 64), xmm3);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+128), xmm2);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+192), xmm7);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+256), xmm1);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+320), xmm11);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+384), xmm6);
    _mm_store_si128((__m128i*)(dprofile+16*j+1536+448), xmm15);
  }

  //  dprofile_dump7(dprofile);
}


// Register usage
// xmm0:  H0
// xmm1:  H1
// xmm2:  H2
// xmm3:  H3
// xmm4:  F0
// xmm5:  F1
// xmm6:  F2
// xmm7:  F3
// xmm8:  N0
// xmm9:  N1
// xmm10: N2
// xmm11: N3
// xmm12: E
// xmm13: temporary
// xmm14: Q 
// xmm15: R


/* 
   initialize Q, R, H0-3, F0-3, 
   loop index (r11), 
   query length (loop end double) (r10), loop end single (r12) 
*/


#define INITIALIZE				\
  "        movq      %3, rax             \n"	\
  "        movdqa    (rax), xmm14        \n"	\
  "        movq      %4, rax             \n"	\
  "        movdqa    (rax), xmm15        \n"	\
  "        movq      %10, rax            \n"	\
  "        movdqa    (rax), xmm0         \n"	\
  "        movdqa    (%7), xmm7          \n"    \
  "        movdqa    xmm7, xmm4          \n"    \
  "        movdqa    xmm7, xmm1          \n"    \
  "        paddusb   xmm15, xmm7         \n"    \
  "        movdqa    xmm7, xmm5          \n"    \
  "        movdqa    xmm7, xmm2          \n"    \
  "        paddusb   xmm15, xmm7         \n"    \
  "        movdqa    xmm7, xmm6          \n"    \
  "        movdqa    xmm7, xmm3          \n"    \
  "        paddusb   xmm15, xmm7         \n"    \
  "        movq      %5, r12             \n"	\
  "        shlq      $3, r12             \n"	\
  "        movq      r12, r10            \n"	\
  "        andq      $-16, r10           \n"	\
  "        xorq      r11, r11            \n" 



#define ONESTEP(H, N, F, V, UP, LEFT)		\
  "        paddusb   "V", "H"            \n"	\
  "        pminub    "F", "H"            \n"	\
  "        pminub    xmm12, "H"          \n"	\
  "        movdqa    "H", xmm13          \n"    \
  "        pcmpeqb   "F", xmm13          \n"    \
  "        pmovmskb  xmm13, edx          \n"	\
  "        movw      dx, "UP"            \n"	\
  "        movdqa    "H", xmm13          \n"	\
  "        pcmpeqb   xmm12, xmm13        \n"    \
  "        pmovmskb  xmm13, edx          \n"	\
  "        movw      dx, "LEFT"          \n"	\
  "        paddusb   xmm15, "F"          \n"	\
  "        paddusb   xmm15, xmm12        \n"	\
  "        movdqa    "H", "N"            \n"	\
  "        paddusb   xmm14, "H"          \n"	\
  "        pminub    "H", xmm12          \n"	\
  "        pminub    "H", "F"            \n"


inline void donormal8(__m128i * Sm,
		      __m128i * hep,
		      __m128i ** qp,
		      __m128i * Qm,
		      __m128i * Rm,
		      long ql,
		      __m128i * Zm,
		      __m128i * F0,
		      WORD * Up,
		      WORD * Left,
		      __m128i * H0
		      )
{
#ifdef DEBUG
  printf("donormal\n");
  printf("Sm=%p\n", Sm);
  printf("hep=%p\n", hep);
  printf("qp=%p\n", qp);
  printf("Qm=%p\n", Qm);
  printf("Rm=%p\n", Rm);
  printf("qlen=%ld\n", ql);
  printf("Zm=%p\n", Zm);
#endif
  
  __asm__
    __volatile__
    ( 
     ".att_syntax noprefix    # Change assembler syntax \n"
     
     INITIALIZE
     
     "        jmp       2f                  \n"
     
     "1:      movq      0(%2,r11,1), rax    \n" // load x from qp[qi]
     "        movdqa    0(%1,r11,4), xmm8   \n" // load N0
     "        movdqa    16(%1,r11,4), xmm12 \n" // load E
     
     ONESTEP("xmm0", "xmm9",        "xmm4", " 0(rax)", " 0(%8,r11,1)", " 0(%9,r11,1)")
     ONESTEP("xmm1", "xmm10",       "xmm5", "16(rax)", " 2(%8,r11,1)", " 2(%9,r11,1)")
     ONESTEP("xmm2", "xmm11",       "xmm6", "32(rax)", " 4(%8,r11,1)", " 4(%9,r11,1)")
     ONESTEP("xmm3", "0(%1,r11,4)", "xmm7", "48(rax)", " 6(%8,r11,1)", " 6(%9,r11,1)")
     
     "        movdqa    xmm12, 16(%1,r11,4) \n" // save E
     "        movq      8(%2,r11,1), rax    \n" // load x from qp[qi+1]
     "        movdqa    32(%1,r11,4), xmm0  \n" // load H0
     "        movdqa    48(%1,r11,4), xmm12 \n" // load E
     
     ONESTEP("xmm8",  "xmm1",           "xmm4", "0(rax)" , " 8(%8,r11,1)", " 8(%9,r11,1)")
     ONESTEP("xmm9",  "xmm2",           "xmm5", "16(rax)", "10(%8,r11,1)", "10(%9,r11,1)")
     ONESTEP("xmm10", "xmm3",           "xmm6", "32(rax)", "12(%8,r11,1)", "12(%9,r11,1)")
     ONESTEP("xmm11", "32(%1,r11,4)",   "xmm7", "48(rax)", "14(%8,r11,1)", "14(%9,r11,1)")
     
     "        movdqa    xmm12, 48(%1,r11,4) \n" // save E
     "        addq      $16, r11            \n" // qi++
     "2:      cmpq      r11, r10            \n" // qi = ql4 ?
     "        jne       1b                  \n" // loop
     
     "4:      cmpq      r11, r12            \n" 
     "        je        3f                  \n"
     "        movq      0(%2,r11,1), rax    \n" // load x from qp[qi]
     "        movdqa    16(%1,r11,4), xmm12 \n" // load E
     
     ONESTEP("xmm0",  "xmm9",          "xmm4", "0(rax)" , " 0(%8,r11,1)", " 0(%9,r11,1)")
     ONESTEP("xmm1",  "xmm10",         "xmm5", "16(rax)", " 2(%8,r11,1)", " 2(%9,r11,1)")
     ONESTEP("xmm2",  "xmm11",         "xmm6", "32(rax)", " 4(%8,r11,1)", " 4(%9,r11,1)")
     ONESTEP("xmm3",  "0(%1,r11,4)",   "xmm7", "48(rax)", " 6(%8,r11,1)", " 6(%9,r11,1)")
     
     "        movdqa    xmm12, 16(%1,r11,4) \n" // save E
     
     "        movdqa    xmm9, xmm1          \n"
     "        movdqa    xmm10, xmm2         \n"
     "        movdqa    xmm11, xmm3         \n"
     "        movdqa    0(%1,r11,4), xmm4   \n"
     "        jmp       5f                  \n"
     
     "3:      movdqa    -32(%1,r11,4), xmm4 \n"
     
     "5:      movq      %0, rax             \n" // save final Hs
     "        movdqa    xmm1, (rax)         \n"
     "        addq      $16, rax            \n"
     "        movdqa    xmm2, (rax)         \n"
     "        addq      $16, rax            \n"
     "        movdqa    xmm3, (rax)         \n"
     "        addq      $16, rax            \n"
     "        movdqa    xmm4, (rax)         \n"
     
     "        .att_syntax prefix      # Change back to standard syntax"
     
     : 
     : "m"(Sm), "r"(hep),  "r"(qp), "m"(Qm), 
       "m"(Rm), "r"(ql),   "m"(Zm), "r"(F0),
       "r"(Up), "r"(Left), 
       "m"(H0)
       
     : "xmm0",  "xmm1",  "xmm2",  "xmm3",
       "xmm4",  "xmm5",  "xmm6",  "xmm7",
       "xmm8",  "xmm9",  "xmm10", "xmm11", 
       "xmm12", "xmm13", "xmm14", "xmm15",
       "rax",   "r10",   "r11",   "r12",
       "rdx",   "cc"
      );
}

inline void domasked8(__m128i * Sm,
		      __m128i * hep,
		      __m128i ** qp,
		      __m128i * Qm, 
		      __m128i * Rm, 
		      long ql,      
		      __m128i * Zm,
		      __m128i * F0,
		      WORD * Up,
		      WORD * Left,
		      __m128i * H0,
		      __m128i * Mm,
		      __m128i * MQ,
		      __m128i * MR)
{
  
#ifdef DEBUG
  printf("domasked\n");
  printf("Sm=%p\n", Sm);
  printf("hep=%p\n", hep);
  printf("qp=%p\n", qp);
  printf("Qm=%p\n", Qm);
  printf("Rm=%p\n", Rm);
  printf("qlen=%ld\n", ql);
  printf("Zm=%p\n", Zm);
  printf("Mm=%p\n", Mm);
#endif

#if 1
  __asm__
    __volatile__
    (
     ".att_syntax noprefix    # Change assembler syntax \n"
     
     INITIALIZE

     "        jmp       2f                   \n"
     
     "1:      movq      0(%2,r11,1), rax     \n" // load x from qp[qi]
     "        movdqa    0(%1,r11,4), xmm8    \n" // load N0
     "        movdqa    16(%1,r11,4), xmm12  \n" // load E
     "        movdqa    (%12), xmm13         \n" 
     "        psubusb   (%11), xmm8          \n" // mask N0
     "        psubusb   (%11), xmm12         \n" // mask E
     "        paddusb   xmm13, xmm8          \n" // init N0
     "        paddusb   xmm13, xmm12         \n" // init E
     "        paddusb   (%13), xmm13         \n" // update
     "        movdqa    xmm13, (%12)         \n"
     
     ONESTEP("xmm0",  "xmm9",          "xmm4", "0(rax)" , "0(%8,r11,1)", "0(%9,r11,1)")
     ONESTEP("xmm1",  "xmm10",         "xmm5", "16(rax)", "2(%8,r11,1)", "2(%9,r11,1)")
     ONESTEP("xmm2",  "xmm11",         "xmm6", "32(rax)", "4(%8,r11,1)", "4(%9,r11,1)")
     ONESTEP("xmm3",  "0(%1,r11,4)",   "xmm7", "48(rax)", "6(%8,r11,1)", "6(%9,r11,1)")
     
     "        movdqa    xmm12, 16(%1,r11,4)  \n" // save E

     "        movq      8(%2,r11,1), rax     \n" // load x from qp[qi+1]
     "        movdqa    32(%1,r11,4), xmm0   \n" // load H0
     "        movdqa    48(%1,r11,4), xmm12  \n" // load E
     "        movdqa    (%12), xmm13         \n"
     "        psubusb   (%11), xmm0          \n" // mask H0
     "        psubusb   (%11), xmm12         \n" // mask E
     "        paddusb   xmm13, xmm0          \n"
     "        paddusb   xmm13, xmm12         \n"
     "        paddusb   (%13), xmm13         \n"
     "        movdqa    xmm13, (%12)         \n"
     
     ONESTEP("xmm8",  "xmm1",           "xmm4", "0(rax)" , " 8(%8,r11,1)", " 8(%9,r11,1)")
     ONESTEP("xmm9",  "xmm2",           "xmm5", "16(rax)", "10(%8,r11,1)", "10(%9,r11,1)")
     ONESTEP("xmm10", "xmm3",           "xmm6", "32(rax)", "12(%8,r11,1)", "12(%9,r11,1)")
     ONESTEP("xmm11", "32(%1,r11,4)",   "xmm7", "48(rax)", "14(%8,r11,1)", "14(%9,r11,1)")
     
     "        movdqa    xmm12, 48(%1,r11,4)  \n" // save E
     "        addq      $16, r11             \n" // qi++
     "2:      cmpq      r11, r10             \n" // qi = ql4 ?
     "        jne       1b                   \n" // loop
     
     "        cmpq      r11, r12             \n" 
     "        je        3f                   \n"
     "        movq      0(%2,r11,1), rax     \n" // load x from qp[qi]
     "        movdqa    16(%1,r11,4), xmm12  \n" // load E
     "        movdqa    (%12), xmm13         \n"
     "        psubusb   (%11), xmm12         \n" // mask E
     "        paddusb   xmm13, xmm12         \n"
     "        paddusb   (%13), xmm13         \n"
     "        movdqa    xmm13, (%12)         \n"
     
     ONESTEP("xmm0",  "xmm9",          "xmm4", "0(rax)" , "0(%8,r11,1)", "0(%9,r11,1)")
     ONESTEP("xmm1",  "xmm10",         "xmm5", "16(rax)", "2(%8,r11,1)", "2(%9,r11,1)")
     ONESTEP("xmm2",  "xmm11",         "xmm6", "32(rax)", "4(%8,r11,1)", "4(%9,r11,1)")
     ONESTEP("xmm3",  "0(%1,r11,4)",   "xmm7", "48(rax)", "6(%8,r11,1)", "6(%9,r11,1)")
     
     "        movdqa    xmm12, 16(%1,r11,4)  \n" // save E
     
     "        movdqa    xmm9, xmm1           \n"
     "        movdqa    xmm10, xmm2          \n"
     "        movdqa    xmm11, xmm3          \n"
     "        movdqa    0(%1,r11,4), xmm4    \n"
     "        jmp       5f                   \n"
     
     "3:      movdqa    -32(%1,r11,4), xmm4  \n"
     
     "5:      movq      %0, rax              \n" // save final Hs
     "        movdqa    xmm1, (rax)          \n"
     "        addq      $16, rax             \n"
     "        movdqa    xmm2, (rax)          \n"
     "        addq      $16, rax             \n"
     "        movdqa    xmm3, (rax)          \n"
     "        addq      $16, rax             \n"
     "        movdqa    xmm4, (rax)          \n"
     
     "        .att_syntax prefix      # Change back to standard syntax"
     
     : 
     
     : "m"(Sm), "r"(hep),"r"(qp), "m"(Qm), 
       "m"(Rm), "r"(ql), "m"(Zm), "r"(F0),
       "r"(Up), "r"(Left),
       "m"(H0), "r"(Mm), "r"(MQ), "r"(MR)
       
     : "xmm0",  "xmm1",  "xmm2",  "xmm3",
       "xmm4",  "xmm5",  "xmm6",  "xmm7",
       "xmm8",  "xmm9",  "xmm10", "xmm11", 
       "xmm12", "xmm13", "xmm14", "xmm15",
       "rax",   "r10",   "r11",   "r12",
       "rdx",   "cc"
     );
#endif
}

unsigned long backtrack(char * qseq,
			char * dseq,
			unsigned long qlen,
			unsigned long dlen,
			WORD * up_array, 
			WORD * left_array,
			unsigned long offset,
			unsigned long size,
			int channel)
{
  unsigned long diff = 0;

#if 0

  printf("Dumping backtracking array\n");

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long index = (offset + longestdbsequence*4*(j/4) + 4*i + (j&3)) % size;
      unsigned long mask = 1 << channel;

      //      printf("c, o, s, i, j, index=%d, %lu, %lu, %lu, %lu, %ld\n", channel, offset, size, i, j, index);

      if (up_array[index] & mask)
      {
	if (left_array[index] & mask)
	  printf("+");
	else
	  printf("^");
      }
      else if (left_array[index] & mask)
      {
	printf("<");
      }
      else
      {
	printf("\\");
      }
    }
    printf("\n");
  }
  
#endif

  //#define SHOWALN

  long i = qlen - 1;
  long j = dlen - 1;

  while ((i>=0) && (j>=0))
  {
    unsigned long index = (offset + longestdbsequence*4*(j/4) + 4*i + (j&3)) % size;
    unsigned long mask = 1 << channel;

    if (up_array[index] & mask)
    {
      diff++;
      i--;
#ifdef SHOWALN
      printf("D");
#endif
    }
    else if (left_array[index] & mask)
    {
      diff++;
      j--;
#ifdef SHOWALN
      printf("I");
#endif
    }
    else
    {
      if (qseq[i] == dseq[j])
      {
#ifdef SHOWALN
        printf("=");
#endif
      }
      else
      {
#ifdef SHOWALN
        printf("X");
#endif
        diff++;
      }
      i--;
      j--;
    }
  }

#ifdef SHOWALN

  while(i>=0)
  {
    diff++;
    i--;
    printf("D");
  }

  while(j>=0)
  {
    diff++;
    j--;
    printf("I");
  }

  printf("\n");

#else

  diff += i + j + 2;

#endif

  return diff;
}

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
	     WORD * up_array,
	     WORD * left_array)
{
  __m128i Q, R, T, M, T0, MQ, MR;
  __m128i *hep, **qp;

  BYTE * d_begin[CHANNELS];
  BYTE * d_end[CHANNELS];
  unsigned long d_offset[CHANNELS];
  BYTE * d_address[CHANNELS];
  unsigned long d_length[CHANNELS];
  
  __m128i dseqalloc[CDEPTH];
  
  __m128i H0, F0;
  __m128i S[4];

  BYTE * dseq = (BYTE*) & dseqalloc;
  BYTE zero;

  long seq_id[CHANNELS];
  unsigned long next_id = 0;
  unsigned long done;
  
  //memset(hearray, 0x00, qlen * 32);

  T0 = _mm_set_epi8(0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 
		    0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff);
  Q  = _mm_set1_epi8(gap_open_penalty+gap_extend_penalty);
  R  = _mm_set1_epi8(gap_extend_penalty);

  zero = 0;
  done = 0;

  hep = (__m128i*) hearray;
  qp = (__m128i**) q_start;

#ifdef DEBUG
  //  printf("Searching %ld sequences...\n", sequences);
#endif

  for (int c=0; c<CHANNELS; c++)
  {
    d_begin[c] = &zero;
    d_end[c] = d_begin[c];
    seq_id[c] = -1;
    /*
    d_offset[c] = 0;
    d_length[c] = 0;
    d_address[c] = 0;
    */
  }
  
  /*
  S[0] = _mm_setzero_si128();
  S[1] = _mm_setzero_si128();
  S[2] = _mm_setzero_si128();
  S[3] = _mm_setzero_si128();
  */

  F0 = _mm_setzero_si128();
  H0 = _mm_setzero_si128();
  
  int easy = 0;

  WORD * up = up_array;
  WORD * left = left_array;

  //  printf("dirbs, qlen, longest: %ld, %ld, %lu\n", dirbufferlen, qlen, longestdbsequence);

  while(1)
  {

    if (easy)
    {
      // fill all channels

      for(int c=0; c<CHANNELS; c++)
      {
	for(int j=0; j<CDEPTH; j++)
	{
	  if (d_begin[c] < d_end[c])
	    dseq[CHANNELS*j+c] = *(d_begin[c]++);
	  else
	    dseq[CHANNELS*j+c] = 0;
	}
	if (d_begin[c] == d_end[c])
	  easy = 0;
      }

      dprofile_shuffle8(dprofile, score_matrix, dseq);
      
      donormal8(S, hep, qp, &Q, &R, qlen, 0, &F0, up, left, &H0);
    }
    else
    {
      // One or more sequences ended in the previous block 
      // We have to switch over to a new sequence

      easy = 1;

      M = _mm_setzero_si128();
      T = T0;
      for (int c=0; c<CHANNELS; c++)
      {
	if (d_begin[c] < d_end[c])
	{
	  // this channel has more sequence

	  for(int j=0; j<CDEPTH; j++)
	  {
	    if (d_begin[c] < d_end[c])
	      dseq[CHANNELS*j+c] = *(d_begin[c]++);
	    else
	      dseq[CHANNELS*j+c] = 0;
	  }
	  if (d_begin[c] == d_end[c])
	    easy = 0;
	}
	else
	{
	  // sequence in channel c ended
	  // change of sequence

	  M = _mm_xor_si128(M, T);

	  long cand_id = seq_id[c];
	  
	  
	  if (cand_id >= 0)
	  {
	    // printf("Completed channel %d, sequence %ld\n", c, cand_id);
	    // save score

	    char * dbseq = (char*) d_address[c];
	    long dbseqlen = d_length[c];
	    long z = (dbseqlen+3) % 4;
	    long score = ((BYTE*)S)[z*16+c];
	    scores[cand_id] = score;
	    
	    unsigned long diff;

	    if (score < 255)
	    {
	      long offset = d_offset[c];
	      diff = backtrack(query.seq, dbseq, qlen, dbseqlen,
			       up_array, left_array,
			       offset,
			       dirbuffersize, c);
	    }
	    else
	    {
	      diff = 255;
	    }

	    diffs[cand_id] = diff;

	    done++;
	  }

	  if (next_id < sequences)
	  {
	    // get next sequence
	    seq_id[c] = next_id;
	    long seqno = seqnos[next_id];
	    char* address;
	    long length;

	    db_getsequenceandlength(seqno, & address, & length);
		      
	    // printf("Seqno: %ld Address: %p\n", seqno, address);
	    d_address[c] = (BYTE*) address;
	    d_length[c] = length;

	    d_begin[c] = (unsigned char*) address;
	    d_end[c] = (unsigned char*) address + length;
	    d_offset[c] = up - up_array;
	    next_id++;
	    
	    ((BYTE*)&H0)[c] = 0;
	    ((BYTE*)&F0)[c] = gap_open_penalty + gap_extend_penalty;
	    
	    
	    // fill channel
	    for(int j=0; j<CDEPTH; j++)
	    {
	      if (d_begin[c] < d_end[c])
		dseq[CHANNELS*j+c] = *(d_begin[c]++);
	      else
		dseq[CHANNELS*j+c] = 0;
	    }
	    if (d_begin[c] == d_end[c])
	      easy = 0;
	  }
	  else
	  {
	    // no more sequences, empty channel
	    seq_id[c] = -1;
	    d_begin[c] = &zero;
	    d_end[c] = d_begin[c];
	    for (int j=0; j<CDEPTH; j++)
	      dseq[CHANNELS*j+c] = 0;
	  }


	}

	T = _mm_slli_si128(T, 1);
      }

      if (done == sequences)
	break;
	  
      dprofile_shuffle8(dprofile, score_matrix, dseq);
	  
      MQ = _mm_and_si128(M, Q); //  // _mm_add_epi8(Q,R));
      MR = _mm_and_si128(M, R);
      
      domasked8(S, hep, qp, &Q, &R, qlen, 0, &F0, up, left, &H0, &M, &MQ, &MR);
    }
    
#if 0
    /* dump H */

    printf("col: %d\n", col);
    col += CDEPTH;

    printf("row: H-values, E-values\n");
    for(int i=0; i<qlen; i++)
    {
      printf("%3d:", i);
      for(int j=0; j<16; j++)
	printf(" %02x", hearray[32*i + j]);
      printf(" -");
      for(int j=0; j<16; j++)
	printf(" %02x", hearray[32*i + j + 16]);
      printf("\n");
    }
#endif

    F0 = _mm_adds_epu8(F0, R);
    F0 = _mm_adds_epu8(F0, R);
    F0 = _mm_adds_epu8(F0, R);
    H0 = F0;
    F0 = _mm_adds_epu8(F0, R);


    /*
    for(int i=0; i<CHANNELS; i++)
    {
      ((BYTE*)&H0)[i] = MIN(((BYTE*)&F0)[i] + 3*gap_extend_penalty, 255);
      ((BYTE*)&F0)[i] = MIN(((BYTE*)&F0)[i] + 4*gap_extend_penalty, 255);
    }
    */

    //    printf("up_array, up: %p %p\n", up_array, up);

    up += 4*longestdbsequence;
    left += 4*longestdbsequence;
    
    if (up >= up_array + dirbuffersize)
    {
      up -= dirbuffersize;
      left -= dirbuffersize;
    }
  }
}
