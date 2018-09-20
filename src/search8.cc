/*
    Swarm

    Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe

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

#define CHANNELS 16
#define CDEPTH 4

#define MATRIXWIDTH 16

void dprofile_dump8(BYTE * dprofile)
{
  char * ss = sym_nt;

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

  //  dprofile_dump8(dprofile);
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

/*
   Sorry for the assembler code below. This code was originally written
   several years ago when compilers were not that good at compiling
   intrinsics to optimal code.
   Similar code using intrinsics instead of assembler is available in
   the vsearch codebase.
*/

#define INITIALIZE                                      \
  "        movq      %3, %%rax               \n"        \
  "        movdqa    (%%rax), %%xmm14        \n"        \
  "        movq      %4, %%rax               \n"        \
  "        movdqa    (%%rax), %%xmm15        \n"        \
  "        movq      %9, %%rax               \n"        \
  "        movdqa    (%%rax), %%xmm0         \n"        \
  "        movdqa    (%7), %%xmm7            \n"        \
  "        movdqa    %%xmm7, %%xmm3          \n"        \
  "        psubusb   %%xmm14, %%xmm3         \n"        \
  "        movdqa    %%xmm3, %%xmm1          \n"        \
  "        paddusb   %%xmm15, %%xmm3         \n"        \
  "        movdqa    %%xmm3, %%xmm2          \n"        \
  "        paddusb   %%xmm15, %%xmm3         \n"        \
  "        movdqa    %%xmm7, %%xmm4          \n"        \
  "        paddusb   %%xmm15, %%xmm7         \n"        \
  "        movdqa    %%xmm7, %%xmm5          \n"        \
  "        paddusb   %%xmm15, %%xmm7         \n"        \
  "        movdqa    %%xmm7, %%xmm6          \n"        \
  "        paddusb   %%xmm15, %%xmm7         \n"        \
  "        movq      %5, %%r12               \n"        \
  "        shlq      $3, %%r12               \n"        \
  "        movq      %%r12, %%r10            \n"        \
  "        andq      $-16, %%r10             \n"        \
  "        xorq      %%r11, %%r11            \n"

#define ONESTEP(H, N, F, V, DIR)                \
  "        paddusb   " V ", " H "          \n"  \
  "        movdqa    " H ", %%xmm13        \n"  \
  "        pminub    " F ", " H "          \n"  \
  "        pcmpeqb   " H ", %%xmm13        \n"  \
  "        pmovmskb  %%xmm13, %%edx        \n"  \
  "        movw      %%dx, 0+" DIR "       \n"  \
  "        pminub    %%xmm12, " H "        \n"  \
  "        movdqa    " H ", %%xmm13        \n"  \
  "        pcmpeqb   %%xmm12, %%xmm13      \n"  \
  "        pmovmskb  %%xmm13, %%edx        \n"  \
  "        movw      %%dx, 2+" DIR "       \n"  \
  "        movdqa    " H ", " N "          \n"  \
  "        paddusb   %%xmm14, " H "        \n"  \
  "        paddusb   %%xmm15, " F "        \n"  \
  "        paddusb   %%xmm15, %%xmm12      \n"  \
  "        movdqa    " H ", %%xmm13        \n"  \
  "        pminub    " H ", " F "          \n"  \
  "        pcmpeqb   " F ", %%xmm13        \n"  \
  "        pmovmskb  %%xmm13, %%edx        \n"  \
  "        movw      %%dx, 4+" DIR "       \n"  \
  "        movdqa    " H ", %%xmm13        \n"  \
  "        pminub    " H ", %%xmm12        \n"  \
  "        pcmpeqb   %%xmm12, %%xmm13      \n"  \
  "        pmovmskb  %%xmm13, %%edx        \n"  \
  "        movw      %%dx, 6+" DIR "       \n"


inline void donormal8(__m128i * Sm,
                      __m128i * hep,
                      __m128i ** qp,
                      __m128i * Qm,
                      __m128i * Rm,
                      long ql,
                      __m128i * Zm,
                      __m128i * F0,
                      unsigned long * dir,
                      __m128i * H0
                      )
{
  __asm__
    __volatile__
    (
     INITIALIZE

     "        jmp       2f                  \n"

     "1:      movq      0(%2,%%r11,1), %%rax    \n" // load x from qp[qi]
     "        movdqa    0(%1,%%r11,4), %%xmm8   \n" // load N0
     "        movdqa    16(%1,%%r11,4), %%xmm12 \n" // load E

     ONESTEP("%%xmm0",  "%%xmm9",         "%%xmm4", " 0(%%rax)", " 0(%8,%%r11,4)")
     ONESTEP("%%xmm1",  "%%xmm10",        "%%xmm5", "16(%%rax)", " 8(%8,%%r11,4)")
     ONESTEP("%%xmm2",  "%%xmm11",        "%%xmm6", "32(%%rax)", "16(%8,%%r11,4)")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",  "%%xmm7", "48(%%rax)", "24(%8,%%r11,4)")

     "        movdqa    %%xmm12, 16(%1,%%r11,4) \n" // save E
     "        movq      8(%2,%%r11,1), %%rax    \n" // load x from qp[qi+1]
     "        movdqa    32(%1,%%r11,4), %%xmm0  \n" // load H0
     "        movdqa    48(%1,%%r11,4), %%xmm12 \n" // load E

     ONESTEP("%%xmm8",  "%%xmm1",         "%%xmm4", " 0(%%rax)", "32(%8,%%r11,4)")
     ONESTEP("%%xmm9",  "%%xmm2",         "%%xmm5", "16(%%rax)", "40(%8,%%r11,4)")
     ONESTEP("%%xmm10", "%%xmm3",         "%%xmm6", "32(%%rax)", "48(%8,%%r11,4)")
     ONESTEP("%%xmm11", "32(%1,%%r11,4)", "%%xmm7", "48(%%rax)", "56(%8,%%r11,4)")

     "        movdqa    %%xmm12, 48(%1,%%r11,4) \n" // save E
     "        addq      $16, %%r11              \n" // qi++
     "2:      cmpq      %%r11, %%r10            \n" // qi = ql4 ?
     "        jne       1b                      \n" // loop

     "4:      cmpq      %%r11, %%r12            \n"
     "        je        3f                      \n"
     "        movq      0(%2,%%r11,1), %%rax    \n" // load x from qp[qi]
     "        movdqa    16(%1,%%r11,4), %%xmm12 \n" // load E

     ONESTEP("%%xmm0",  "%%xmm9",         "%%xmm4", " 0(%%rax)", " 0(%8,%%r11,4)")
     ONESTEP("%%xmm1",  "%%xmm10",        "%%xmm5", "16(%%rax)", " 8(%8,%%r11,4)")
     ONESTEP("%%xmm2",  "%%xmm11",        "%%xmm6", "32(%%rax)", "16(%8,%%r11,4)")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",  "%%xmm7", "48(%%rax)", "24(%8,%%r11,4)")

     "        movdqa    %%xmm12, 16(%1,%%r11,4) \n" // save E

     "        movdqa    %%xmm9, %%xmm1          \n"
     "        movdqa    %%xmm10, %%xmm2         \n"
     "        movdqa    %%xmm11, %%xmm3         \n"
     "        movdqa    0(%1,%%r11,4), %%xmm4   \n"
     "        jmp       5f                      \n"

     "3:      movdqa    -32(%1,%%r11,4), %%xmm4 \n"

     "5:      movq      %0, %%rax               \n" // save final Hs
     "        movdqa    %%xmm1, (%%rax)         \n"
     "        addq      $16, %%rax              \n"
     "        movdqa    %%xmm2, (%%rax)         \n"
     "        addq      $16, %%rax              \n"
     "        movdqa    %%xmm3, (%%rax)         \n"
     "        addq      $16, %%rax              \n"
     "        movdqa    %%xmm4, (%%rax)         \n"

     :
     : "m"(Sm), "r"(hep),  "r"(qp), "m"(Qm),
       "m"(Rm), "r"(ql),   "m"(Zm), "r"(F0),
       "r"(dir),"m"(H0)

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
                      unsigned long * dir,
                      __m128i * H0,
                      __m128i * Mm,
                      __m128i * MQ,
                      __m128i * MR,
                      __m128i * MQ0)
{
  __asm__
    __volatile__
    (
     INITIALIZE

     "        jmp       2f                       \n"

     "1:      movq      0(%2,%%r11,1), %%rax     \n" // load x from qp[qi]
     "        movdqa    0(%1,%%r11,4), %%xmm8    \n" // load N0
     "        movdqa    16(%1,%%r11,4), %%xmm12  \n" // load E
     "        movdqa    (%11), %%xmm13           \n"
     "        psubusb   (%10), %%xmm8            \n" // mask N0
     "        psubusb   (%10), %%xmm12           \n" // mask E
     "        paddusb   %%xmm13, %%xmm8          \n" // init N0
     "        paddusb   %%xmm13, %%xmm12         \n" // init E
     "        paddusb   (%13), %%xmm12           \n" // fix E
     "        paddusb   (%12), %%xmm13           \n" // update
     "        movdqa    %%xmm13, (%11)           \n"

     ONESTEP("%%xmm0",  "%%xmm9",         "%%xmm4", " 0(%%rax)", " 0(%8,%%r11,4)")
     ONESTEP("%%xmm1",  "%%xmm10",        "%%xmm5", "16(%%rax)", " 8(%8,%%r11,4)")
     ONESTEP("%%xmm2",  "%%xmm11",        "%%xmm6", "32(%%rax)", "16(%8,%%r11,4)")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",  "%%xmm7", "48(%%rax)", "24(%8,%%r11,4)")

     "        movdqa    %%xmm12, 16(%1,%%r11,4)  \n" // save E

     "        movq      8(%2,%%r11,1), %%rax     \n" // load x from qp[qi+1]
     "        movdqa    32(%1,%%r11,4), %%xmm0   \n" // load H0
     "        movdqa    48(%1,%%r11,4), %%xmm12  \n" // load E
     "        movdqa    (%11), %%xmm13           \n"
     "        psubusb   (%10), %%xmm0            \n" // mask H0
     "        psubusb   (%10), %%xmm12           \n" // mask E
     "        paddusb   %%xmm13, %%xmm0          \n"
     "        paddusb   %%xmm13, %%xmm12         \n"
     "        paddusb   (%13), %%xmm12           \n" // fix E
     "        paddusb   (%12), %%xmm13           \n"
     "        movdqa    %%xmm13, (%11)           \n"

     ONESTEP("%%xmm8",  "%%xmm1",         "%%xmm4", " 0(%%rax)", "32(%8,%%r11,4)")
     ONESTEP("%%xmm9",  "%%xmm2",         "%%xmm5", "16(%%rax)", "40(%8,%%r11,4)")
     ONESTEP("%%xmm10", "%%xmm3",         "%%xmm6", "32(%%rax)", "48(%8,%%r11,4)")
     ONESTEP("%%xmm11", "32(%1,%%r11,4)", "%%xmm7", "48(%%rax)", "56(%8,%%r11,4)")

     "        movdqa    %%xmm12, 48(%1,%%r11,4)  \n" // save E
     "        addq      $16, %%r11               \n" // qi++
     "2:      cmpq      %%r11, %%r10             \n" // qi = ql4 ?
     "        jne       1b                       \n" // loop

     "        cmpq      %%r11, %%r12             \n"
     "        je        3f                       \n"
     "        movq      0(%2,%%r11,1), %%rax     \n" // load x from qp[qi]
     "        movdqa    16(%1,%%r11,4), %%xmm12  \n" // load E
     "        movdqa    (%11), %%xmm13           \n"
     "        psubusb   (%10), %%xmm12           \n" // mask E
     "        paddusb   %%xmm13, %%xmm12         \n"
     "        paddusb   (%13), %%xmm12           \n" // fix E
     "        paddusb   (%12), %%xmm13           \n"
     "        movdqa    %%xmm13, (%11)           \n"

     ONESTEP("%%xmm0",  "%%xmm9",          "%%xmm4", "0(%%rax)" , " 0(%8,%%r11,4)")
     ONESTEP("%%xmm1",  "%%xmm10",         "%%xmm5", "16(%%rax)", " 8(%8,%%r11,4)")
     ONESTEP("%%xmm2",  "%%xmm11",         "%%xmm6", "32(%%rax)", "16(%8,%%r11,4)")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",   "%%xmm7", "48(%%rax)", "24(%8,%%r11,4)")

     "        movdqa    %%xmm12, 16(%1,%%r11,4)  \n" // save E

     "        movdqa    %%xmm9, %%xmm1           \n"
     "        movdqa    %%xmm10, %%xmm2          \n"
     "        movdqa    %%xmm11, %%xmm3          \n"
     "        movdqa    0(%1,%%r11,4), %%xmm4    \n"
     "        jmp       5f                       \n"

     "3:      movdqa    -32(%1,%%r11,4), %%xmm4  \n"

     "5:      movq      %0, %%rax                \n" // save final Hs
     "        movdqa    %%xmm1, (%%rax)          \n"
     "        addq      $16, %%rax               \n"
     "        movdqa    %%xmm2, (%%rax)          \n"
     "        addq      $16, %%rax               \n"
     "        movdqa    %%xmm3, (%%rax)          \n"
     "        addq      $16, %%rax               \n"
     "        movdqa    %%xmm4, (%%rax)          \n"

     :

     : "m"(Sm), "r"(hep),"r"(qp), "m"(Qm),
       "m"(Rm), "r"(ql), "m"(Zm), "r"(F0),
       "r"(dir),"m"(H0),
       "r"(Mm), "r"(MQ), "r"(MR), "r"(MQ0)

     : "xmm0",  "xmm1",  "xmm2",  "xmm3",
       "xmm4",  "xmm5",  "xmm6",  "xmm7",
       "xmm8",  "xmm9",  "xmm10", "xmm11",
       "xmm12", "xmm13", "xmm14", "xmm15",
       "rax",   "r10",   "r11",   "r12",
       "rdx",   "cc"
     );
}

unsigned long backtrack(char * qseq,
                        char * dseq,
                        unsigned long qlen,
                        unsigned long dlen,
                        unsigned long * dirbuffer,
                        unsigned long offset,
                        unsigned long dirbuffersize,
                        unsigned long channel,
                        unsigned long * alignmentlengthp)
{
  unsigned long maskup      = 1UL << (channel+ 0);
  unsigned long maskleft    = 1UL << (channel+16);
  unsigned long maskextup   = 1UL << (channel+32);
  unsigned long maskextleft = 1UL << (channel+48);

#if 0

  printf("Dumping backtracking array\n");

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4) + 4*i + (j&3)) % dirbuffersize];
      if (d & maskleft)
      {
        printf("<");
      }
      else if (!(d & maskup))
      {
          printf("^");
      }
      else
      {
        printf("\\");
      }
    }
    printf("\n");
  }

  printf("Dumping gap extension array\n");

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4)
                                   + 4*i + (j&3)) % dirbuffersize];
      if (!(d & maskextup))
      {
        if (!(d & maskextleft))
          printf("+");
        else
          printf("^");
      }
      else if (!(d & maskextleft))
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

  long i = qlen - 1;
  long j = dlen - 1;
  unsigned long aligned = 0;
  unsigned long matches = 0;
  char op = 0;

#undef SHOWALIGNMENT
#ifdef SHOWALIGNMENT
  printf("alignment, reversed: ");
#endif

  while ((i>=0) && (j>=0))
  {
    aligned++;

    unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4)
                                 + 4*i + (j&3)) % dirbuffersize];

    if ((op == 'I') && (!(d & maskextleft)))
    {
      j--;
    }
    else if ((op == 'D') && (!(d & maskextup)))
    {
      i--;
    }
    else if (d & maskleft)
    {
      j--;
      op = 'I';
    }
    else if (!(d & maskup))
    {
      i--;
      op = 'D';
    }
    else
    {
      if (nt_extract(qseq, i) == nt_extract(dseq, j))
        matches++;
      i--;
      j--;
      op = 'M';
    }

#ifdef SHOWALIGNMENT
    printf("%c", op);
#endif

  }

  while (i>=0)
    {
      aligned++;
      i--;
#ifdef SHOWALIGNMENT
      printf("D");
#endif
    }

  while (j>=0)
    {
      aligned++;
      j--;
#ifdef SHOWALIGNMENT
      printf("I");
#endif
    }

#ifdef SHOWALIGNMENT
  printf("\n");
#endif

  * alignmentlengthp = aligned;
  return aligned - matches;
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
             unsigned long * alignmentlengths,
             unsigned long qlen,
             unsigned long dirbuffersize,
             unsigned long * dirbuffer)
{
  __m128i Q, R, T, M, T0, MQ, MR, MQ0;
  __m128i *hep, **qp;

  unsigned long d_pos[CHANNELS];
  unsigned long d_offset[CHANNELS];
  char * d_address[CHANNELS];
  unsigned long d_length[CHANNELS];

  __m128i dseqalloc[CDEPTH];

  __m128i H0, F0;
  __m128i S[4];

  BYTE * dseq = (BYTE*) & dseqalloc;

  long seq_id[CHANNELS];
  unsigned long next_id = 0;
  unsigned long done;

  T0 = _mm_set_epi8(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1);
  Q  = _mm_set1_epi8(gap_open_penalty+gap_extend_penalty);
  R  = _mm_set1_epi8(gap_extend_penalty);

  done = 0;

  hep = (__m128i*) hearray;
  qp = (__m128i**) q_start;

  for (int c=0; c<CHANNELS; c++)
  {
    d_address[c] = 0;
    d_pos[c] = 0;
    d_length[c] = 0;
    seq_id[c] = -1;
  }

  F0 = _mm_setzero_si128();
  H0 = _mm_setzero_si128();

  int easy = 0;

  unsigned long * dir = dirbuffer;

  while(1)
  {

    if (easy)
    {
      // fill all channels

      for(int c=0; c<CHANNELS; c++)
      {
        for(int j=0; j<CDEPTH; j++)
        {
          if (d_pos[c] < d_length[c])
            dseq[CHANNELS*j+c] = 1 + nt_extract(d_address[c], d_pos[c]++);
          else
            dseq[CHANNELS*j+c] = 0;
        }
        if (d_pos[c] == d_length[c])
          easy = 0;
      }

      if (ssse3_present)
        dprofile_shuffle8(dprofile, score_matrix, dseq);
      else
        dprofile_fill8(dprofile, score_matrix, dseq);

      donormal8(S, hep, qp, &Q, &R, qlen, 0, &F0, dir, &H0);
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
        if (d_pos[c] < d_length[c])
        {
          // this channel has more sequence

          for(int j=0; j<CDEPTH; j++)
          {
            if (d_pos[c] < d_length[c])
              dseq[CHANNELS*j+c] = 1 + nt_extract(d_address[c], d_pos[c]++);
            else
              dseq[CHANNELS*j+c] = 0;
          }
          if (d_pos[c] == d_length[c])
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
                               dirbuffer,
                               offset,
                               dirbuffersize, c,
                               alignmentlengths + cand_id);
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

            d_address[c] = address;
            d_length[c] = length;

            d_pos[c] = 0;
            d_offset[c] = dir - dirbuffer;
            next_id++;

            ((BYTE*)&H0)[c] = 0;
            ((BYTE*)&F0)[c] = 2 * gap_open_penalty + 2 * gap_extend_penalty;

            // fill channel
            for(int j=0; j<CDEPTH; j++)
            {
              if (d_pos[c] < d_length[c])
                dseq[CHANNELS*j+c] = 1 + nt_extract(d_address[c], d_pos[c]++);
              else
                dseq[CHANNELS*j+c] = 0;
            }
            if (d_pos[c] == d_length[c])
              easy = 0;
          }
          else
          {
            // no more sequences, empty channel
            seq_id[c] = -1;
            d_address[c] = 0;
            d_pos[c] = 0;
            d_length[c] = 0;
            for (int j=0; j<CDEPTH; j++)
              dseq[CHANNELS*j+c] = 0;
          }


        }

        T = _mm_slli_si128(T, 1);
      }

      if (done == sequences)
        break;

      if (ssse3_present)
        dprofile_shuffle8(dprofile, score_matrix, dseq);
      else
        dprofile_fill8(dprofile, score_matrix, dseq);

      MQ = _mm_and_si128(M, Q);
      MR = _mm_and_si128(M, R);
      MQ0 = MQ;

      domasked8(S, hep, qp, &Q, &R, qlen, 0, &F0, dir, &H0, &M, &MQ, &MR,
                &MQ0);
    }

    F0 = _mm_adds_epu8(F0, R);
    F0 = _mm_adds_epu8(F0, R);
    F0 = _mm_adds_epu8(F0, R);
    H0 = _mm_subs_epu8(F0, Q);
    F0 = _mm_adds_epu8(F0, R);


    dir += 4*longestdbsequence;
    if (dir >= dirbuffer + dirbuffersize)
      dir -= dirbuffersize;
  }
}
