/*
    SWARM

    Copyright (C) 2012-2014 Torbjorn Rognes and Frederic Mahe

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

#define CHANNELS 8
#define CDEPTH 4

#define SHUFFLE 1

void dprofile_dump16(WORD * dprofile)
{
  char * s = sym_nt;
  printf("\ndprofile:\n");
  for(int i=0; i<32; i++)
  {
    printf("%c: ",s[i]);
    for(int k=0; k<CDEPTH; k++)
    {
      printf("[");
      for(int j=0; j<CHANNELS; j++)
	printf(" %3d", (short) dprofile[CHANNELS*CDEPTH*i + CHANNELS*k + j]);
      printf("]");
    }
    printf("\n");
  }
  exit(1);
}

inline void dprofile_shuffle16(WORD * dprofile,
			       WORD * score_matrix,
			       BYTE * dseq_byte)
{
  __m128i m0, m1, m2, m3, t0, t1, t2, t3, t4, t5, u0, u1, u2, u3, u4, zero, one;
  __m128i * dseq = (__m128i*) dseq_byte;

  // instructions? 32 + 5*12 = 92

  // Requires SSSE3 due to pshufb

  zero = _mm_setzero_si128();
  one  = _mm_set1_epi16(1);
  
  t0 = _mm_load_si128(dseq+0);
  
  m0 = _mm_unpacklo_epi8(t0, zero);
  m0 = _mm_slli_epi16(m0, 1);
  t1 = _mm_adds_epu16(m0, one);
  t1 = _mm_slli_epi16(t1, 8);
  m0 = _mm_or_si128(m0, t1);
  
  m1 = _mm_unpackhi_epi8(t0, zero);
  m1 = _mm_slli_epi16(m1, 1);
  t2 = _mm_adds_epu16(m1, one);
  t2 = _mm_slli_epi16(t2, 8);
  m1 = _mm_or_si128(m1, t2);

  t3 = _mm_load_si128(dseq+1);
  
  m2 = _mm_unpacklo_epi8(t3, zero);
  m2 = _mm_slli_epi16(m2, 1);
  t4 = _mm_adds_epu16(m2, one);
  t4 = _mm_slli_epi16(t4, 8);
  m2 = _mm_or_si128(m2, t4);

  m3 = _mm_unpackhi_epi8(t3, zero);
  m3 = _mm_slli_epi16(m3, 1);
  t5 = _mm_adds_epu16(m3, one);
  t5 = _mm_slli_epi16(t5, 8);
  m3 = _mm_or_si128(m3, t5);
  
#define profline(j)                                     \
  u0 = _mm_load_si128((__m128i*)(score_matrix)+4*j);    \
  u1 = _mm_shuffle_epi8(u0, m0);                        \
  u2 = _mm_shuffle_epi8(u0, m1);                        \
  u3 = _mm_shuffle_epi8(u0, m2);                        \
  u4 = _mm_shuffle_epi8(u0, m3);                        \
  _mm_store_si128((__m128i*)(dprofile)+4*j+0, u1);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+1, u2);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+2, u3);      \
  _mm_store_si128((__m128i*)(dprofile)+4*j+3, u4)

  profline(0);
  profline(1);
  profline(2);
  profline(3);
  profline(4);

#if 0
  dprofile_dump16(dprofile);
#endif
}



inline void dprofile_fill16(WORD * dprofile_word,
			    WORD * score_matrix_word,
			    BYTE * dseq)
{
  __m128i xmm0,  xmm1,  xmm2,  xmm3,  xmm4,  xmm5,  xmm6,  xmm7;
  __m128i xmm8,  xmm9,  xmm10, xmm11, xmm12, xmm13, xmm14, xmm15;
  __m128i xmm16, xmm17, xmm18, xmm19, xmm20, xmm21, xmm22, xmm23;
  __m128i xmm24, xmm25, xmm26, xmm27, xmm28, xmm29, xmm30, xmm31;

  // clocks? 4*(8+3*(8+4)+8) = 52*4 = 208
  
  for (int j=0; j<CDEPTH; j++)
  {
    int d[CHANNELS];
    for(int z=0; z<CHANNELS; z++)
      d[z] = dseq[j*CHANNELS+z] << 5;
      
    // for(int i=0; i<24; i += 8)
    // for(int i=0; i<32; i += 8)
    for(int i=0; i<8; i += 8)
    {
      xmm0  = _mm_load_si128((__m128i*)(score_matrix_word + d[0] + i));
      xmm1  = _mm_load_si128((__m128i*)(score_matrix_word + d[1] + i));
      xmm2  = _mm_load_si128((__m128i*)(score_matrix_word + d[2] + i));
      xmm3  = _mm_load_si128((__m128i*)(score_matrix_word + d[3] + i));
      xmm4  = _mm_load_si128((__m128i*)(score_matrix_word + d[4] + i));
      xmm5  = _mm_load_si128((__m128i*)(score_matrix_word + d[5] + i));
      xmm6  = _mm_load_si128((__m128i*)(score_matrix_word + d[6] + i));
      xmm7  = _mm_load_si128((__m128i*)(score_matrix_word + d[7] + i));
      
      xmm8  = _mm_unpacklo_epi16(xmm0,  xmm1);
      xmm9  = _mm_unpackhi_epi16(xmm0,  xmm1);
      xmm10 = _mm_unpacklo_epi16(xmm2,  xmm3);
      xmm11 = _mm_unpackhi_epi16(xmm2,  xmm3);
      xmm12 = _mm_unpacklo_epi16(xmm4,  xmm5);
      xmm13 = _mm_unpackhi_epi16(xmm4,  xmm5);
      xmm14 = _mm_unpacklo_epi16(xmm6,  xmm7);
      xmm15 = _mm_unpackhi_epi16(xmm6,  xmm7);
      
      xmm16 = _mm_unpacklo_epi32(xmm8,  xmm10);
      xmm17 = _mm_unpackhi_epi32(xmm8,  xmm10);
      xmm18 = _mm_unpacklo_epi32(xmm12, xmm14);
      xmm19 = _mm_unpackhi_epi32(xmm12, xmm14);
      xmm20 = _mm_unpacklo_epi32(xmm9,  xmm11);
      xmm21 = _mm_unpackhi_epi32(xmm9,  xmm11);
      xmm22 = _mm_unpacklo_epi32(xmm13, xmm15);
      xmm23 = _mm_unpackhi_epi32(xmm13, xmm15);
      
      xmm24 = _mm_unpacklo_epi64(xmm16, xmm18);
      xmm25 = _mm_unpackhi_epi64(xmm16, xmm18);
      xmm26 = _mm_unpacklo_epi64(xmm17, xmm19);
      xmm27 = _mm_unpackhi_epi64(xmm17, xmm19);
      xmm28 = _mm_unpacklo_epi64(xmm20, xmm22);
      xmm29 = _mm_unpackhi_epi64(xmm20, xmm22);
      xmm30 = _mm_unpacklo_epi64(xmm21, xmm23);
      xmm31 = _mm_unpackhi_epi64(xmm21, xmm23);
      
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+0) + CHANNELS*j), xmm24);
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+1) + CHANNELS*j), xmm25);
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+2) + CHANNELS*j), xmm26);
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+3) + CHANNELS*j), xmm27);
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+4) + CHANNELS*j), xmm28);
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+5) + CHANNELS*j), xmm29);
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+6) + CHANNELS*j), xmm30);
      _mm_store_si128((__m128i*)(dprofile_word + CDEPTH*CHANNELS*(i+7) + CHANNELS*j), xmm31);
    }
  }
#if 0
  dprofile_dump16(dprofile_word);
#endif
}

// Due to the use of pminsw instead of pminuw below, the code works
// only with 15-bit values

#define INITIALIZE					\
  "        movq      %3, %%rax               \n"	\
  "        movdqa    (%%rax), %%xmm14        \n"	\
  "        movq      %4, %%rax               \n"	\
  "        movdqa    (%%rax), %%xmm15        \n"	\
  "        movq      %9, %%rax               \n"	\
  "        movdqa    (%%rax), %%xmm0         \n"	\
  "        movdqa    (%7), %%xmm7            \n"	\
  "        movdqa    %%xmm7, %%xmm4          \n"	\
  "        movdqa    %%xmm7, %%xmm1          \n"	\
  "        paddusw   %%xmm15, %%xmm7         \n"	\
  "        movdqa    %%xmm7, %%xmm5          \n"	\
  "        movdqa    %%xmm7, %%xmm2          \n"	\
  "        paddusw   %%xmm15, %%xmm7         \n"	\
  "        movdqa    %%xmm7, %%xmm6          \n"	\
  "        movdqa    %%xmm7, %%xmm3          \n"	\
  "        paddusw   %%xmm15, %%xmm7         \n"	\
  "        movq      %5, %%r12               \n"	\
  "        shlq      $3, %%r12               \n"	\
  "        movq      %%r12, %%r10            \n"	\
  "        andq      $-16, %%r10             \n"	\
  "        xorq      %%r11, %%r11            \n" 

#define ONESTEP(H, N, F, V, DIR)			\
  "        paddusw   "V", "H"                \n"	\
  "        movdqa    "H", %%xmm13            \n"	\
  "        pcmpgtw   "F", %%xmm13            \n"	\
  "        pmovmskb  %%xmm13, %%edx          \n"	\
  "        movw      %%dx, 0+"DIR"           \n"	\
  "        pminsw    "F", "H"                \n"	\
  "        pminsw    %%xmm12, "H"            \n"	\
  "        movdqa    "H", %%xmm13            \n"	\
  "        pcmpeqw   %%xmm12, %%xmm13        \n"	\
  "        pmovmskb  %%xmm13, %%edx          \n"	\
  "        movw      %%dx, 2+"DIR"           \n"	\
  "        movdqa    "H", "N"                \n"	\
  "        paddusw   %%xmm14, "H"            \n"	\
  "        paddusw   %%xmm15, "F"            \n"	\
  "        paddusw   %%xmm15, %%xmm12        \n"	\
  "        movdqa    "H", %%xmm13            \n"	\
  "        pcmpgtw   "F", %%xmm13            \n"	\
  "        pmovmskb  %%xmm13, %%edx          \n"	\
  "        movw      %%dx, 4+"DIR"           \n"	\
  "        movdqa    "H", %%xmm13            \n"	\
  "        pcmpgtw   %%xmm12, %%xmm13        \n"	\
  "        pmovmskb  %%xmm13, %%edx          \n"	\
  "        movw      %%dx, 6+"DIR"           \n"	\
  "        pminsw    "H", %%xmm12            \n"	\
  "        pminsw    "H", "F"                \n"


inline void donormal16(__m128i * Sm,
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
     
     ONESTEP("%%xmm0", "%%xmm9",        "%%xmm4", " 0(%%rax)", " 0(%8,%%r11,4)")
     ONESTEP("%%xmm1", "%%xmm10",       "%%xmm5", "16(%%rax)", " 8(%8,%%r11,4)")
     ONESTEP("%%xmm2", "%%xmm11",       "%%xmm6", "32(%%rax)", "16(%8,%%r11,4)")
     ONESTEP("%%xmm3", "0(%1,%%r11,4)", "%%xmm7", "48(%%rax)", "24(%8,%%r11,4)")
     
     "        movdqa    %%xmm12, 16(%1,%%r11,4) \n" // save E
     "        movq      8(%2,%%r11,1), %%rax    \n" // load x from qp[qi+1]
     "        movdqa    32(%1,%%r11,4), %%xmm0  \n" // load H0
     "        movdqa    48(%1,%%r11,4), %%xmm12 \n" // load E
     
     ONESTEP("%%xmm8",  "%%xmm1",           "%%xmm4", "0(%%rax)" , "32(%8,%%r11,4)")
     ONESTEP("%%xmm9",  "%%xmm2",           "%%xmm5", "16(%%rax)", "40(%8,%%r11,4)")
     ONESTEP("%%xmm10", "%%xmm3",           "%%xmm6", "32(%%rax)", "48(%8,%%r11,4)")
     ONESTEP("%%xmm11", "32(%1,%%r11,4)",   "%%xmm7", "48(%%rax)", "56(%8,%%r11,4)")
     
     "        movdqa    %%xmm12, 48(%1,%%r11,4) \n" // save E
     "        addq      $16, %%r11              \n" // qi++
     "2:      cmpq      %%r11, %%r10            \n" // qi = ql4 ?
     "        jne       1b                      \n" // loop
     
     "4:      cmpq      %%r11, %%r12            \n" 
     "        je        3f                      \n"
     "        movq      0(%2,%%r11,1), %%rax    \n" // load x from qp[qi]
     "        movdqa    16(%1,%%r11,4), %%xmm12 \n" // load E
     
     ONESTEP("%%xmm0",  "%%xmm9",          "%%xmm4", "0(%%rax)" , " 0(%8,%%r11,4)")
     ONESTEP("%%xmm1",  "%%xmm10",         "%%xmm5", "16(%%rax)", " 8(%8,%%r11,4)")
     ONESTEP("%%xmm2",  "%%xmm11",         "%%xmm6", "32(%%rax)", "16(%8,%%r11,4)")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",   "%%xmm7", "48(%%rax)", "24(%8,%%r11,4)")
     
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

inline void domasked16(__m128i * Sm,
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
		       __m128i * MR)
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
     "        psubusw   (%10), %%xmm8            \n" // mask N0
     "        psubusw   (%10), %%xmm12           \n" // mask E
     "        paddusw   %%xmm13, %%xmm8          \n" // init N0
     "        paddusw   %%xmm13, %%xmm12         \n" // init E
     "        paddusw   (%12), %%xmm13           \n" // update
     "        movdqa    %%xmm13, (%11)           \n"
     
     ONESTEP("%%xmm0",  "%%xmm9",          "%%xmm4", "0(%%rax)" , " 0(%8,%%r11,4)")
     ONESTEP("%%xmm1",  "%%xmm10",         "%%xmm5", "16(%%rax)", " 8(%8,%%r11,4)")
     ONESTEP("%%xmm2",  "%%xmm11",         "%%xmm6", "32(%%rax)", "16(%8,%%r11,4)")
     ONESTEP("%%xmm3",  "0(%1,%%r11,4)",   "%%xmm7", "48(%%rax)", "24(%8,%%r11,4)")
     
     "        movdqa    %%xmm12, 16(%1,%%r11,4)  \n" // save E

     "        movq      8(%2,%%r11,1), %%rax     \n" // load x from qp[qi+1]
     "        movdqa    32(%1,%%r11,4), %%xmm0   \n" // load H0
     "        movdqa    48(%1,%%r11,4), %%xmm12  \n" // load E
     "        movdqa    (%11), %%xmm13           \n"
     "        psubusw   (%10), %%xmm0            \n" // mask H0
     "        psubusw   (%10), %%xmm12           \n" // mask E
     "        paddusw   %%xmm13, %%xmm0          \n"
     "        paddusw   %%xmm13, %%xmm12         \n"
     "        paddusw   (%12), %%xmm13           \n"
     "        movdqa    %%xmm13, (%11)           \n"
     
     ONESTEP("%%xmm8",  "%%xmm1",           "%%xmm4", "0(%%rax)" , "32(%8,%%r11,4)")
     ONESTEP("%%xmm9",  "%%xmm2",           "%%xmm5", "16(%%rax)", "40(%8,%%r11,4)")
     ONESTEP("%%xmm10", "%%xmm3",           "%%xmm6", "32(%%rax)", "48(%8,%%r11,4)")
     ONESTEP("%%xmm11", "32(%1,%%r11,4)",   "%%xmm7", "48(%%rax)", "56(%8,%%r11,4)")
     
     "        movdqa    %%xmm12, 48(%1,%%r11,4)  \n" // save E
     "        addq      $16, %%r11               \n" // qi++
     "2:      cmpq      %%r11, %%r10             \n" // qi = ql4 ?
     "        jne       1b                       \n" // loop
     
     "        cmpq      %%r11, %%r12             \n" 
     "        je        3f                       \n"
     "        movq      0(%2,%%r11,1), %%rax     \n" // load x from qp[qi]
     "        movdqa    16(%1,%%r11,4), %%xmm12  \n" // load E
     "        movdqa    (%11), %%xmm13           \n"
     "        psubusw   (%10), %%xmm12           \n" // mask E
     "        paddusw   %%xmm13, %%xmm12         \n"
     "        paddusw   (%12), %%xmm13           \n"
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
       "r"(dir),
       "m"(H0), "r"(Mm), "r"(MQ), "r"(MR)
       
     : "xmm0",  "xmm1",  "xmm2",  "xmm3",
       "xmm4",  "xmm5",  "xmm6",  "xmm7",
       "xmm8",  "xmm9",  "xmm10", "xmm11", 
       "xmm12", "xmm13", "xmm14", "xmm15",
       "rax",   "r10",   "r11",   "r12",
       "rdx",   "cc"
     );
}

unsigned long backtrack16(char * qseq,
			  char * dseq,
			  unsigned long qlen,
			  unsigned long dlen,
			  unsigned long * dirbuffer,
			  unsigned long offset,
			  unsigned long dirbuffersize,
			  unsigned long channel,
			  unsigned long * alignmentlengthp)
{
  unsigned long maskup      = 3UL << (2*channel+ 0);
  unsigned long maskleft    = 3UL << (2*channel+16);
  unsigned long maskextup   = 3UL << (2*channel+32);
  unsigned long maskextleft = 3UL << (2*channel+48);

#if 0

  printf("Dumping backtracking array\n");

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4) + 4*i + (j&3)) % dir
				  buffersize];
      if (d & maskup)
      {
	if (d & maskleft)
	  printf("+");
	else
	  printf("^");
      }
      else if (d & maskleft)
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

  printf("Dumping gap extension array\n");

  for(unsigned long i=0; i<qlen; i++)
  {
    for(unsigned long j=0; j<dlen; j++)
    {
      unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4) + 4*i + (j&3)) % dir
				  buffersize];
      if (d & maskextup)
      {
	if (d & maskextleft)
	  printf("+");
	else
	  printf("^");
      }
      else if (d & maskextleft)
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

  while ((i>=0) && (j>=0))
  {
    aligned++;

    unsigned long d = 
      dirbuffer[(offset + longestdbsequence*4*(j/4) + 4*i + (j&3)) % dirbuffersize];

    if ((op == 'I') && (d & maskextleft))
    {
      j--;
    }
    else if ((op == 'D') && (d & maskextup))
    {
      i--;
    }
    else if (d & maskleft)
    {
      j--;
      op = 'I';
    }
    else if (d & maskup)
    {
      i--;
      op = 'D';
    }
    else
    {
      if (qseq[i] == dseq[j])
	matches++;
      i--;
      j--;
      op = 'M';
    }
  }
  aligned += i + j + 2;
  * alignmentlengthp = aligned;
  return aligned - matches;
}

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
	      unsigned long * dirbuffer)
{
  __m128i Q, R, T, M, T0, MQ, MR;
  __m128i *hep, **qp;

  BYTE * d_begin[CHANNELS];
  BYTE * d_end[CHANNELS];
  unsigned long d_offset[CHANNELS];
  BYTE * d_address[CHANNELS];
  unsigned long d_length[CHANNELS];
  
  __m128i dseqalloc[CDEPTH];
  
  __m128i H0;
  __m128i F0;
  __m128i S[4];

  BYTE * dseq = (BYTE*) & dseqalloc;
  BYTE zero;

  long seq_id[CHANNELS];
  unsigned long next_id = 0;
  unsigned long done;
  
  T0 = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, 0xffff);
  Q  = _mm_set1_epi16(gap_open_penalty+gap_extend_penalty);
  R  = _mm_set1_epi16(gap_extend_penalty);

  zero = 0;
  done = 0;

  hep = (__m128i*) hearray;
  qp = (__m128i**) q_start;

  for (int c=0; c<CHANNELS; c++)
  {
    d_begin[c] = &zero;
    d_end[c] = d_begin[c];
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
	  if (d_begin[c] < d_end[c])
	    dseq[CHANNELS*j+c] = *(d_begin[c]++);
	  else
	    dseq[CHANNELS*j+c] = 0;
	}
	if (d_begin[c] == d_end[c])
	  easy = 0;
      }

      if (ssse3_present)
	dprofile_shuffle16(dprofile, score_matrix, dseq);
      else
	dprofile_fill16(dprofile, score_matrix, dseq);
      
      donormal16(S, hep, qp, &Q, &R, qlen, 0, &F0, dir, &H0);
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
	    long score = ((WORD*)S)[z*CHANNELS+c];
	    scores[cand_id] = score;
	    
	    unsigned long diff;

	    if (score < 65535)
	      {
		long offset = d_offset[c];
		diff = backtrack16(query.seq, dbseq, qlen, dbseqlen,
				   dirbuffer,
				   offset,
				   dirbuffersize, c,
				   alignmentlengths + cand_id);
	      }
	    else
	      {
		diff = MIN((65535 / penalty_mismatch),
			   (65535 - penalty_gapopen) / penalty_gapextend);
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
	    d_offset[c] = dir - dirbuffer;
	    next_id++;
	    
	    ((WORD*)&H0)[c] = 0;
	    ((WORD*)&F0)[c] = gap_open_penalty + gap_extend_penalty;
	    
	    
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

	T = _mm_slli_si128(T, 2);
      }

      if (done == sequences)
	break;
	  
      if (ssse3_present)
	dprofile_shuffle16(dprofile, score_matrix, dseq);
      else
	dprofile_fill16(dprofile, score_matrix, dseq);
	  
      MQ = _mm_and_si128(M, Q);
      MR = _mm_and_si128(M, R);
      
      domasked16(S, hep, qp, &Q, &R, qlen, 0, &F0, dir, &H0, &M, &MQ, &MR);
    }
    
    F0 = _mm_adds_epu16(F0, R);
    F0 = _mm_adds_epu16(F0, R);
    F0 = _mm_adds_epu16(F0, R);
    H0 = F0;
    F0 = _mm_adds_epu16(F0, R);

    dir += 4*longestdbsequence;
    
    if (dir >= dirbuffer + dirbuffersize)
    {
      dir -= dirbuffersize;
    }
  }
}
