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

#include "swarm.h"

#define CHANNELS 16
#define CDEPTH 4

#ifdef __PPC__

#error Architecture ppcle64 not implemented yet

#elif __aarch64__

#error Architecture aarch64 not implemented yet

#elif __x86_64__

typedef __m128i VECTORTYPE;

#define v_init(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p) _mm_set_epi8(p,o,n,m,l,k,j,i,h,g,f,e,d,c,b,a)
#define v_load(a) _mm_load_si128((VECTORTYPE *)(a))
#define v_load_64(a) _mm_loadl_epi64((VECTORTYPE *)(a))
#define v_store(a, b) _mm_store_si128((VECTORTYPE *)(a), (b))
#define v_merge_lo_8(a, b) _mm_unpacklo_epi8((a),(b))
#define v_merge_hi_8(a, b) _mm_unpackhi_epi8((a),(b))
#define v_merge_lo_16(a, b) _mm_unpacklo_epi16((a),(b))
#define v_merge_hi_16(a, b) _mm_unpackhi_epi16((a),(b))
#define v_merge_lo_32(a, b) _mm_unpacklo_epi32((a),(b))
#define v_merge_hi_32(a, b) _mm_unpackhi_epi32((a),(b))
#define v_merge_lo_64(a, b) _mm_unpacklo_epi64((a),(b))
#define v_merge_hi_64(a, b) _mm_unpackhi_epi64((a),(b))
#define v_min(a, b) _mm_min_epi8((a), (b))
#define v_min_u(a, b) _mm_min_epu8((a), (b))
#define v_add(a, b) _mm_adds_epu8((a), (b))
#define v_sub(a, b) _mm_subs_epu8((a), (b))
#define v_dup(a) _mm_set1_epi8(a)
#define v_zero v_dup(0)
#define v_and(a, b) _mm_and_si128((a), (b))
#define v_xor(a, b) _mm_xor_si128((a), (b))
#define v_shift_left(a) _mm_slli_si128((a), 1)
#define v_mask_gt(a, b) _mm_movemask_epi8(_mm_cmpgt_epi8((a), (b)))
#define v_mask_eq(a, b) _mm_movemask_epi8(_mm_cmpeq_epi8((a), (b)))

#else

#error Unknown Architecture

#endif

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
  __m128i reg0,  reg1, reg2,  reg3,  reg4,  reg5,  reg6,  reg7;
  __m128i reg8,  reg9, reg10, reg11, reg12, reg13, reg14, reg15;

  for(int j=0; j<CDEPTH; j++)
    {
      unsigned d[CHANNELS];
      for(int i=0; i<CHANNELS; i++)
        d[i] = dseq[j*CHANNELS+i] << 5;

      reg0  = v_load_64(score_matrix + d[ 0]);
      reg2  = v_load_64(score_matrix + d[ 2]);
      reg4  = v_load_64(score_matrix + d[ 4]);
      reg6  = v_load_64(score_matrix + d[ 6]);
      reg8  = v_load_64(score_matrix + d[ 8]);
      reg10 = v_load_64(score_matrix + d[10]);
      reg12 = v_load_64(score_matrix + d[12]);
      reg14 = v_load_64(score_matrix + d[14]);

      reg0  = v_merge_lo_8(reg0,  *(VECTORTYPE*)(score_matrix + d[ 1]));
      reg2  = v_merge_lo_8(reg2,  *(VECTORTYPE*)(score_matrix + d[ 3]));
      reg4  = v_merge_lo_8(reg4,  *(VECTORTYPE*)(score_matrix + d[ 5]));
      reg6  = v_merge_lo_8(reg6,  *(VECTORTYPE*)(score_matrix + d[ 7]));
      reg8  = v_merge_lo_8(reg8,  *(VECTORTYPE*)(score_matrix + d[ 9]));
      reg10 = v_merge_lo_8(reg10, *(VECTORTYPE*)(score_matrix + d[11]));
      reg12 = v_merge_lo_8(reg12, *(VECTORTYPE*)(score_matrix + d[13]));
      reg14 = v_merge_lo_8(reg14, *(VECTORTYPE*)(score_matrix + d[15]));

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile+16*j+  0, reg0);
      v_store(dprofile+16*j+ 64, reg3);
      v_store(dprofile+16*j+128, reg2);
      v_store(dprofile+16*j+192, reg7);
      v_store(dprofile+16*j+256, reg1);
      v_store(dprofile+16*j+320, reg11);
      v_store(dprofile+16*j+384, reg6);
      v_store(dprofile+16*j+448, reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(score_matrix + 8 + d[0 ]);
      reg1  = v_load_64(score_matrix + 8 + d[1 ]);
      reg2  = v_load_64(score_matrix + 8 + d[2 ]);
      reg3  = v_load_64(score_matrix + 8 + d[3 ]);
      reg4  = v_load_64(score_matrix + 8 + d[4 ]);
      reg5  = v_load_64(score_matrix + 8 + d[5 ]);
      reg6  = v_load_64(score_matrix + 8 + d[6 ]);
      reg7  = v_load_64(score_matrix + 8 + d[7 ]);
      reg8  = v_load_64(score_matrix + 8 + d[8 ]);
      reg9  = v_load_64(score_matrix + 8 + d[9 ]);
      reg10 = v_load_64(score_matrix + 8 + d[10]);
      reg11 = v_load_64(score_matrix + 8 + d[11]);
      reg12 = v_load_64(score_matrix + 8 + d[12]);
      reg13 = v_load_64(score_matrix + 8 + d[13]);
      reg14 = v_load_64(score_matrix + 8 + d[14]);
      reg15 = v_load_64(score_matrix + 8 + d[15]);

      reg0  = v_merge_lo_8(reg0,  reg1);
      reg2  = v_merge_lo_8(reg2,  reg3);
      reg4  = v_merge_lo_8(reg4,  reg5);
      reg6  = v_merge_lo_8(reg6,  reg7);
      reg8  = v_merge_lo_8(reg8,  reg9);
      reg10 = v_merge_lo_8(reg10, reg11);
      reg12 = v_merge_lo_8(reg12, reg13);
      reg14 = v_merge_lo_8(reg14, reg15);

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile+16*j+512+  0, reg0);
      v_store(dprofile+16*j+512+ 64, reg3);
      v_store(dprofile+16*j+512+128, reg2);
      v_store(dprofile+16*j+512+192, reg7);
      v_store(dprofile+16*j+512+256, reg1);
      v_store(dprofile+16*j+512+320, reg11);
      v_store(dprofile+16*j+512+384, reg6);
      v_store(dprofile+16*j+512+448, reg15);


      reg0  = v_load_64(score_matrix + 16 + d[0 ]);
      reg2  = v_load_64(score_matrix + 16 + d[2 ]);
      reg4  = v_load_64(score_matrix + 16 + d[4 ]);
      reg6  = v_load_64(score_matrix + 16 + d[6 ]);
      reg8  = v_load_64(score_matrix + 16 + d[8 ]);
      reg10 = v_load_64(score_matrix + 16 + d[10]);
      reg12 = v_load_64(score_matrix + 16 + d[12]);
      reg14 = v_load_64(score_matrix + 16 + d[14]);

      reg0  = v_merge_lo_8(reg0,  *(VECTORTYPE*)(score_matrix + 16 + d[ 1]));
      reg2  = v_merge_lo_8(reg2,  *(VECTORTYPE*)(score_matrix + 16 + d[ 3]));
      reg4  = v_merge_lo_8(reg4,  *(VECTORTYPE*)(score_matrix + 16 + d[ 5]));
      reg6  = v_merge_lo_8(reg6,  *(VECTORTYPE*)(score_matrix + 16 + d[ 7]));
      reg8  = v_merge_lo_8(reg8,  *(VECTORTYPE*)(score_matrix + 16 + d[ 9]));
      reg10 = v_merge_lo_8(reg10, *(VECTORTYPE*)(score_matrix + 16 + d[11]));
      reg12 = v_merge_lo_8(reg12, *(VECTORTYPE*)(score_matrix + 16 + d[13]));
      reg14 = v_merge_lo_8(reg14, *(VECTORTYPE*)(score_matrix + 16 + d[15]));

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile+16*j+1024+  0, reg0);
      v_store(dprofile+16*j+1024+ 64, reg3);
      v_store(dprofile+16*j+1024+128, reg2);
      v_store(dprofile+16*j+1024+192, reg7);
      v_store(dprofile+16*j+1024+256, reg1);
      v_store(dprofile+16*j+1024+320, reg11);
      v_store(dprofile+16*j+1024+384, reg6);
      v_store(dprofile+16*j+1024+448, reg15);


      // loads not aligned on 16 byte boundary, cannot load and unpack in one instr.

      reg0  = v_load_64(score_matrix + 24 + d[ 0]);
      reg1  = v_load_64(score_matrix + 24 + d[ 1]);
      reg2  = v_load_64(score_matrix + 24 + d[ 2]);
      reg3  = v_load_64(score_matrix + 24 + d[ 3]);
      reg4  = v_load_64(score_matrix + 24 + d[ 4]);
      reg5  = v_load_64(score_matrix + 24 + d[ 5]);
      reg6  = v_load_64(score_matrix + 24 + d[ 6]);
      reg7  = v_load_64(score_matrix + 24 + d[ 7]);
      reg8  = v_load_64(score_matrix + 24 + d[ 8]);
      reg9  = v_load_64(score_matrix + 24 + d[ 9]);
      reg10 = v_load_64(score_matrix + 24 + d[10]);
      reg11 = v_load_64(score_matrix + 24 + d[11]);
      reg12 = v_load_64(score_matrix + 24 + d[12]);
      reg13 = v_load_64(score_matrix + 24 + d[13]);
      reg14 = v_load_64(score_matrix + 24 + d[14]);
      reg15 = v_load_64(score_matrix + 24 + d[15]);

      reg0  = v_merge_lo_8(reg0,  reg1);
      reg2  = v_merge_lo_8(reg2,  reg3);
      reg4  = v_merge_lo_8(reg4,  reg5);
      reg6  = v_merge_lo_8(reg6,  reg7);
      reg8  = v_merge_lo_8(reg8,  reg9);
      reg10 = v_merge_lo_8(reg10, reg11);
      reg12 = v_merge_lo_8(reg12, reg13);
      reg14 = v_merge_lo_8(reg14, reg15);

      reg1 = reg0;
      reg0 = v_merge_lo_16(reg0, reg2);
      reg1 = v_merge_hi_16(reg1, reg2);
      reg5 = reg4;
      reg4 = v_merge_lo_16(reg4, reg6);
      reg5 = v_merge_hi_16(reg5, reg6);
      reg9 = reg8;
      reg8 = v_merge_lo_16(reg8, reg10);
      reg9 = v_merge_hi_16(reg9, reg10);
      reg13 = reg12;
      reg12 = v_merge_lo_16(reg12, reg14);
      reg13 = v_merge_hi_16(reg13, reg14);

      reg2  = reg0;
      reg0  = v_merge_lo_32(reg0, reg4);
      reg2  = v_merge_hi_32(reg2, reg4);
      reg6  = reg1;
      reg1  = v_merge_lo_32(reg1, reg5);
      reg6  = v_merge_hi_32(reg6, reg5);
      reg10 = reg8;
      reg8  = v_merge_lo_32(reg8, reg12);
      reg10 = v_merge_hi_32(reg10, reg12);
      reg14 = reg9;
      reg9  = v_merge_lo_32(reg9, reg13);
      reg14 = v_merge_hi_32(reg14, reg13);

      reg3  = reg0;
      reg0  = v_merge_lo_64(reg0, reg8);
      reg3  = v_merge_hi_64(reg3, reg8);
      reg7  = reg2;
      reg2  = v_merge_lo_64(reg2, reg10);
      reg7  = v_merge_hi_64(reg7, reg10);
      reg11 = reg1;
      reg1  = v_merge_lo_64(reg1, reg9);
      reg11 = v_merge_hi_64(reg11, reg9);
      reg15 = reg6;
      reg6  = v_merge_lo_64(reg6, reg14);
      reg15 = v_merge_hi_64(reg15, reg14);

      v_store(dprofile+16*j+1536+  0, reg0);
      v_store(dprofile+16*j+1536+ 64, reg3);
      v_store(dprofile+16*j+1536+128, reg2);
      v_store(dprofile+16*j+1536+192, reg7);
      v_store(dprofile+16*j+1536+256, reg1);
      v_store(dprofile+16*j+1536+320, reg11);
      v_store(dprofile+16*j+1536+384, reg6);
      v_store(dprofile+16*j+1536+448, reg15);
    }

  //  dprofile_dump8(dprofile);
}

#define ONESTEP(H, N, F, V, DIR, E, QR, R, W)   \
  H = v_add(H, V);                              \
  W = H;                                        \
  H = v_min_u(H, F);                            \
  *((DIR) + 0) = v_mask_eq(W, H);               \
  H = v_min_u(H, E);                            \
  *((DIR) + 1) = v_mask_eq(H, E);               \
  N = H;                                        \
  H = v_add(H, QR);                             \
  F = v_add(F, R);                              \
  E = v_add(E, R);                              \
  F = v_min_u(H, F);                            \
  *((DIR) + 2) = v_mask_eq(H, F);               \
  E = v_min_u(H, E);                            \
  *((DIR) + 3) = v_mask_eq(H, E);

void align_cells_regular_8(VECTORTYPE * Sm,
                           VECTORTYPE * hep,
                           VECTORTYPE ** qp,
                           VECTORTYPE * Qm,
                           VECTORTYPE * Rm,
                           long ql,
                           VECTORTYPE * F0,
                           unsigned long * dir_long,
                           VECTORTYPE * H0)
{
  VECTORTYPE Q, R, E, W;
  VECTORTYPE h0, h1, h2, h3, h4, h5, h6, h7, h8;
  VECTORTYPE f0, f1, f2, f3;
  VECTORTYPE * x;

  unsigned short * dir = (unsigned short *) dir_long;

  Q = *Qm;
  R = *Rm;

  f0 = *F0;
  f1 = v_add(f0, R);
  f2 = v_add(f1, R);
  f3 = v_add(f2, R);

  h0 = *H0;
  h1 = v_sub(f0, Q);
  h2 = v_add(h1, R);
  h3 = v_add(h2, R);

  for(long i = 0; i < ql; i++)
    {
      x = qp[i + 0];
      h4 = hep[2*i + 0];
      E  = hep[2*i + 1];
      ONESTEP(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R, W);
      ONESTEP(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R, W);
      ONESTEP(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R, W);
      ONESTEP(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R, W);
      hep[2*i + 0] = h8;
      hep[2*i + 1] = E;
      h0 = h4;
      h1 = h5;
      h2 = h6;
      h3 = h7;
    }

  Sm[0] = h5;
  Sm[1] = h6;
  Sm[2] = h7;
  Sm[3] = h8;
}

inline void align_cells_masked_8(VECTORTYPE * Sm,
                                 VECTORTYPE * hep,
                                 VECTORTYPE ** qp,
                                 VECTORTYPE * Qm,
                                 VECTORTYPE * Rm,
                                 long ql,
                                 VECTORTYPE * F0,
                                 unsigned long * dir_long,
                                 VECTORTYPE * H0,
                                 VECTORTYPE * Mm,
                                 VECTORTYPE * MQ,
                                 VECTORTYPE * MR,
                                 VECTORTYPE * MQ0)
{
  VECTORTYPE Q, R, E, W;
  VECTORTYPE h0, h1, h2, h3, h4, h5, h6, h7, h8;
  VECTORTYPE f0, f1, f2, f3;
  VECTORTYPE * x;
  long z, i;

  unsigned short * dir = (unsigned short *) dir_long;

  Q = *Qm;
  R = *Rm;

  f0 = *F0;
  f1 = v_add(f0, R);
  f2 = v_add(f1, R);
  f3 = v_add(f2, R);

  h0 = *H0;
  h1 = v_sub(f0, Q);
  h2 = v_add(h1, R);
  h3 = v_add(h2, R);
  h4 = v_zero;

  z = ql - (ql & 1);
  i = 0;
  while (i < z)
    {
      h4 = hep[2*i + 0];
      E  = hep[2*i + 1];
      x = qp[i + 0];

      /* mask h4 and E */
      h4 = v_sub(h4, *Mm);
      E  = v_sub(E,  *Mm);

      /* init h4 and E */
      h4 = v_add(h4, *MQ);
      E  = v_add(E,  *MQ);
      E  = v_add(E,  *MQ0);

      /* update MQ */
      *MQ = v_add(*MQ,  *MR);

      ONESTEP(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R, W);
      ONESTEP(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R, W);
      ONESTEP(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R, W);
      ONESTEP(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R, W);
      hep[2*i + 0] = h8;
      hep[2*i + 1] = E;

      h0 = hep[2*i + 2];
      E  = hep[2*i + 3];
      x = qp[i +  1];

      /* mask h0 and E */
      h0 = v_sub(h0, *Mm);
      E  = v_sub(E,  *Mm);

      /* init h0 and E */
      h0 = v_add(h0, *MQ);
      E  = v_add(E,  *MQ);
      E  = v_add(E,  *MQ0);

      /* update MQ */
      *MQ = v_add(*MQ, *MR);

      ONESTEP(h4, h1, f0, x[0], dir + 16*i + 16, E, Q, R, W);
      ONESTEP(h5, h2, f1, x[1], dir + 16*i + 20, E, Q, R, W);
      ONESTEP(h6, h3, f2, x[2], dir + 16*i + 24, E, Q, R, W);
      ONESTEP(h7, h4, f3, x[3], dir + 16*i + 28, E, Q, R, W);
      hep[2*i + 2] = h4;
      hep[2*i + 3] = E;

      i += 2;
    }

  if (i < ql)
    {
      E  = hep[2*i + 1];
      x = qp[i + 0];

      /* mask E */
      E  = v_sub(E,  *Mm);

      /* init E */
      E  = v_add(E,  *MQ);
      E  = v_add(E,  *MQ0);

      /* update MQ */
      *MQ = v_add(*MQ,  *MR);

      ONESTEP(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R, W);
      ONESTEP(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R, W);
      ONESTEP(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R, W);
      ONESTEP(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R, W);
      hep[2*i + 0] = h8;
      hep[2*i + 1] = E;

      Sm[0] = h5;
      Sm[1] = h6;
      Sm[2] = h7;
      Sm[3] = h8;
    }
  else
    {
      Sm[0] = h1;
      Sm[1] = h2;
      Sm[2] = h3;
      Sm[3] = h4;
    }
}

inline unsigned long backtrack_8(char * qseq,
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
          unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4)
                                       + 4*i + (j&3)) % dirbuffersize];
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

#define SHOWALIGNMENT
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
  VECTORTYPE Q, R, T, M, T0, MQ, MR, MQ0;
  VECTORTYPE *hep, **qp;

  unsigned long d_pos[CHANNELS];
  unsigned long d_offset[CHANNELS];
  char * d_address[CHANNELS];
  unsigned long d_length[CHANNELS];

  VECTORTYPE dseqalloc[CDEPTH];

  VECTORTYPE H0;
  VECTORTYPE F0;
  VECTORTYPE S[4];

  BYTE * dseq = (BYTE*) & dseqalloc;

  long seq_id[CHANNELS];
  unsigned long next_id = 0;
  unsigned long done;

  T0 = v_init(-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
  Q  = v_dup(gap_open_penalty+gap_extend_penalty);
  R  = v_dup(gap_extend_penalty);

  done = 0;

  hep = (VECTORTYPE*) hearray;
  qp = (VECTORTYPE**) q_start;

  for (int c=0; c<CHANNELS; c++)
    {
      d_address[c] = 0;
      d_pos[c] = 0;
      d_length[c] = 0;
      seq_id[c] = -1;
    }

  F0 = v_zero;
  H0 = v_zero;

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

          align_cells_regular_8(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
        }
      else
        {
          // One or more sequences ended in the previous block
          // We have to switch over to a new sequence

          easy = 1;

          M = v_zero;
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

                  M = v_xor(M, T);

                  long cand_id = seq_id[c];

                  if (cand_id >= 0)
                    {
                      // save score

                      char * dbseq = (char*) d_address[c];
                      long dbseqlen = d_length[c];
                      long z = (dbseqlen+3) % 4;
                      long score = ((BYTE*)S)[z*CHANNELS+c];
                      scores[cand_id] = score;

                      unsigned long diff;

                      if (score < 255)
                        {
                          long offset = d_offset[c];
                          diff = backtrack_8(query.seq, dbseq, qlen, dbseqlen,
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

              T = v_shift_left(T);
            }

          if (done == sequences)
            break;

          if (ssse3_present)
            dprofile_shuffle8(dprofile, score_matrix, dseq);
          else
            dprofile_fill8(dprofile, score_matrix, dseq);

          MQ = v_and(M, Q);
          MR = v_and(M, R);
          MQ0 = MQ;

          align_cells_masked_8(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
        }

      F0 = v_add(F0, R);
      F0 = v_add(F0, R);
      F0 = v_add(F0, R);
      H0 = v_sub(F0, Q);
      F0 = v_add(F0, R);

      dir += 4*longestdbsequence;
      if (dir >= dirbuffer + dirbuffersize)
        dir -= dirbuffersize;
    }
}
