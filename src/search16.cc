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

#define CHANNELS 8
#define CDEPTH 4

#ifdef __PPC__

#error Architecture ppcle64 not implemented yet

#elif __aarch64__

#error Architecture aarch64 not implemented yet

#elif __x86_64__

typedef signed short VECTORELEMENTTYPE;
typedef __m128i VECTORTYPE;

#define v_init(a,b,c,d,e,f,g,h) _mm_set_epi16(h,g,f,e,d,c,b,a)
#define v_load(a) _mm_load_si128((VECTORTYPE *)(a))
#define v_store(a, b) _mm_store_si128((VECTORTYPE *)(a), (b))
#define v_merge_lo_16(a, b) _mm_unpacklo_epi16((a),(b))
#define v_merge_hi_16(a, b) _mm_unpackhi_epi16((a),(b))
#define v_merge_lo_32(a, b) _mm_unpacklo_epi32((a),(b))
#define v_merge_hi_32(a, b) _mm_unpackhi_epi32((a),(b))
#define v_merge_lo_64(a, b) _mm_unpacklo_epi64((a),(b))
#define v_merge_hi_64(a, b) _mm_unpackhi_epi64((a),(b))
#define v_add(a, b) _mm_adds_epi16((a), (b))
#define v_sub(a, b) _mm_subs_epi16((a), (b))
#define v_max(a, b) _mm_max_epi16((a), (b))
#define v_min(a, b) _mm_min_epi16((a), (b))
#define v_add_u(a, b) _mm_adds_epu16((a), (b))
#define v_sub_u(a, b) _mm_subs_epu16((a), (b))
#define v_max_u(a, b) _mm_max_epu16((a), (b))
#define v_min_u(a, b) _mm_min_epu16((a), (b))
#define v_dup(a) _mm_set1_epi16(a)
#define v_zero v_dup(0)
#define v_and(a, b) _mm_and_si128((a), (b))
#define v_xor(a, b) _mm_xor_si128((a), (b))
#define v_shift_left(a) _mm_slli_si128((a), 2)
#define v_mask_gt(a, b) _mm_movemask_epi8(_mm_cmpgt_epi16((a), (b)))
#define v_mask_eq(a, b) _mm_movemask_epi8(_mm_cmpeq_epi16((a), (b)))

#else

#error Unknown Architecture

#endif

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


inline void dprofile_fill16(WORD * dprofile_word,
                            WORD * score_matrix_word,
                            BYTE * dseq)
{
  VECTORTYPE reg0,  reg1,  reg2,  reg3,  reg4,  reg5,  reg6,  reg7;
  VECTORTYPE reg8,  reg9,  reg10, reg11, reg12, reg13, reg14, reg15;
  VECTORTYPE reg16, reg17, reg18, reg19, reg20, reg21, reg22, reg23;
  VECTORTYPE reg24, reg25, reg26, reg27, reg28, reg29, reg30, reg31;

  for (int j=0; j<CDEPTH; j++)
  {
    int d[CHANNELS];
    for(int z=0; z<CHANNELS; z++)
      d[z] = dseq[j*CHANNELS+z] << 5;

    for(int i=0; i<8; i += 8)
    {
      reg0  = v_load(score_matrix_word + d[0] + i);
      reg1  = v_load(score_matrix_word + d[1] + i);
      reg2  = v_load(score_matrix_word + d[2] + i);
      reg3  = v_load(score_matrix_word + d[3] + i);
      reg4  = v_load(score_matrix_word + d[4] + i);
      reg5  = v_load(score_matrix_word + d[5] + i);
      reg6  = v_load(score_matrix_word + d[6] + i);
      reg7  = v_load(score_matrix_word + d[7] + i);

      reg8  = v_merge_lo_16(reg0,  reg1);
      reg9  = v_merge_hi_16(reg0,  reg1);
      reg10 = v_merge_lo_16(reg2,  reg3);
      reg11 = v_merge_hi_16(reg2,  reg3);
      reg12 = v_merge_lo_16(reg4,  reg5);
      reg13 = v_merge_hi_16(reg4,  reg5);
      reg14 = v_merge_lo_16(reg6,  reg7);
      reg15 = v_merge_hi_16(reg6,  reg7);

      reg16 = v_merge_lo_32(reg8,  reg10);
      reg17 = v_merge_hi_32(reg8,  reg10);
      reg18 = v_merge_lo_32(reg12, reg14);
      reg19 = v_merge_hi_32(reg12, reg14);
      reg20 = v_merge_lo_32(reg9,  reg11);
      reg21 = v_merge_hi_32(reg9,  reg11);
      reg22 = v_merge_lo_32(reg13, reg15);
      reg23 = v_merge_hi_32(reg13, reg15);

      reg24 = v_merge_lo_64(reg16, reg18);
      reg25 = v_merge_hi_64(reg16, reg18);
      reg26 = v_merge_lo_64(reg17, reg19);
      reg27 = v_merge_hi_64(reg17, reg19);
      reg28 = v_merge_lo_64(reg20, reg22);
      reg29 = v_merge_hi_64(reg20, reg22);
      reg30 = v_merge_lo_64(reg21, reg23);
      reg31 = v_merge_hi_64(reg21, reg23);

      v_store(dprofile_word + CDEPTH*CHANNELS*(i+0) + CHANNELS*j, reg24);
      v_store(dprofile_word + CDEPTH*CHANNELS*(i+1) + CHANNELS*j, reg25);
      v_store(dprofile_word + CDEPTH*CHANNELS*(i+2) + CHANNELS*j, reg26);
      v_store(dprofile_word + CDEPTH*CHANNELS*(i+3) + CHANNELS*j, reg27);
      v_store(dprofile_word + CDEPTH*CHANNELS*(i+4) + CHANNELS*j, reg28);
      v_store(dprofile_word + CDEPTH*CHANNELS*(i+5) + CHANNELS*j, reg29);
      v_store(dprofile_word + CDEPTH*CHANNELS*(i+6) + CHANNELS*j, reg30);
      v_store(dprofile_word + CDEPTH*CHANNELS*(i+7) + CHANNELS*j, reg31);
    }
  }
#if 0
  dprofile_dump16(dprofile_word);
#endif
}

/* The code works only with 15-bit values */

#define ONESTEP(H, N, F, V, DIR, E, QR, R)                              \
  H = v_add_u(H, V);                                                    \
  *(DIR+0) = v_mask_gt(H, F);                                           \
  H = v_min(H, F);                                                      \
  H = v_min(H, E);                                                      \
  *(DIR+1) = v_mask_eq(H, E);                                           \
  N = H;                                                                \
  H = v_add_u(H, QR);                                                   \
  F = v_add_u(F, R);                                                    \
  E = v_add_u(E, R);                                                    \
  *(DIR+2) = v_mask_gt(H, F);                                           \
  *(DIR+3) = v_mask_gt(H, E);                                           \
  E = v_min(H, E);                                                      \
  F = v_min(H, F);

inline void donormal16(VECTORTYPE * Sm,
                       VECTORTYPE * hep,
                       VECTORTYPE ** qp,
                       VECTORTYPE * Qm,
                       VECTORTYPE * Rm,
                       long ql,
                       VECTORTYPE * F0,
                       unsigned long * dir_long,
                       VECTORTYPE * H0)
{
  VECTORTYPE Q, R, E;
  VECTORTYPE h0, h1, h2, h3, h4, h5, h6, h7, h8;
  VECTORTYPE f0, f1, f2, f3;
  VECTORTYPE * x;
  long z, i;

  VECTORELEMENTTYPE * dir = (VECTORELEMENTTYPE *) dir_long;

  Q = *Qm;
  R = *Rm;

  f0 = *F0;
  f1 = v_add_u(f0, R);
  f2 = v_add_u(f1, R);
  f3 = v_add_u(f2, R);

  h0 = *H0;
  h1 = v_sub_u(f0, Q);
  h2 = v_add_u(h1, R);
  h3 = v_add_u(h2, R);
  h4 = v_zero;

  z = ql - (ql & 1);
  i = 0;
  while (i < z)
    {
      h4 = hep[2*i + 0];
      E  = hep[2*i + 1];
      x = qp[i + 0];
      ONESTEP(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R);
      ONESTEP(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R);
      ONESTEP(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R);
      ONESTEP(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R);
      hep[2*i + 0] = h8;
      hep[2*i + 1] = E;

      h0 = hep[2*i + 2];
      E  = hep[2*i + 3];
      x = qp[i +  1];
      ONESTEP(h4, h1, f0, x[0], dir + 16*i + 16, E, Q, R);
      ONESTEP(h5, h2, f1, x[1], dir + 16*i + 20, E, Q, R);
      ONESTEP(h6, h3, f2, x[2], dir + 16*i + 24, E, Q, R);
      ONESTEP(h7, h4, f3, x[3], dir + 16*i + 28, E, Q, R);
      hep[2*i + 2] = h4;
      hep[2*i + 3] = E;

      i += 2;
    }

  if (i < ql)
    {
      E  = hep[2*i + 1];
      x = qp[i + 0];
      ONESTEP(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R);
      ONESTEP(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R);
      ONESTEP(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R);
      ONESTEP(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R);
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

inline void domasked16(VECTORTYPE * Sm,
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
  VECTORTYPE Q, R, E;
  VECTORTYPE h0, h1, h2, h3, h4, h5, h6, h7, h8;
  VECTORTYPE f0, f1, f2, f3;
  VECTORTYPE * x;
  long z, i;

  VECTORELEMENTTYPE * dir = (VECTORELEMENTTYPE *) dir_long;

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
      h4 = v_sub_u(h4, *Mm);
      E  = v_sub_u(E,  *Mm);

      /* init h4 and E */
      h4 = v_add_u(h4, *MQ);
      E  = v_add_u(E,  *MQ);
      E  = v_add_u(E,  *MQ0);

      /* update MQ */
      *MQ = v_add_u(*MQ,  *MR);

      ONESTEP(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R);
      ONESTEP(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R);
      ONESTEP(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R);
      ONESTEP(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R);
      hep[2*i + 0] = h8;
      hep[2*i + 1] = E;

      h0 = hep[2*i + 2];
      E  = hep[2*i + 3];
      x = qp[i +  1];

      /* mask h0 and E */
      h0 = v_sub_u(h0, *Mm);
      E  = v_sub_u(E,  *Mm);

      /* init h0 and E */
      h0 = v_add_u(h0, *MQ);
      E  = v_add_u(E,  *MQ);
      E  = v_add_u(E,  *MQ0);

      /* update MQ */
      *MQ = v_add_u(*MQ, *MR);

      ONESTEP(h4, h1, f0, x[0], dir + 16*i + 16, E, Q, R);
      ONESTEP(h5, h2, f1, x[1], dir + 16*i + 20, E, Q, R);
      ONESTEP(h6, h3, f2, x[2], dir + 16*i + 24, E, Q, R);
      ONESTEP(h7, h4, f3, x[3], dir + 16*i + 28, E, Q, R);
      hep[2*i + 2] = h4;
      hep[2*i + 3] = E;

      i += 2;
    }

  if (i < ql)
    {
      E  = hep[2*i + 1];
      x = qp[i + 0];

      /* mask E */
      E  = v_sub_u(E,  *Mm);

      /* init E */
      E  = v_add_u(E,  *MQ);
      E  = v_add_u(E,  *MQ0);

      /* update MQ */
      *MQ = v_add_u(*MQ,  *MR);

      ONESTEP(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R);
      ONESTEP(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R);
      ONESTEP(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R);
      ONESTEP(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R);
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
      unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4)
                                   + 4*i + (j&3)) % dirbuffersize];
      if (d & maskleft)
      {
        printf("<");
      }
      else if (d & maskup)
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

#undef SHOWALIGNMENT
#ifdef SHOWALIGNMENT
  printf("alignment, reversed: ");
#endif

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

  T0 = v_init(-1, 0, 0, 0, 0, 0, 0, 0);
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
        dprofile_shuffle16(dprofile, score_matrix, dseq);
      else
        dprofile_fill16(dprofile, score_matrix, dseq);

      donormal16(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
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

            d_address[c] = address;
            d_length[c] = length;

            d_pos[c] = 0;
            d_offset[c] = dir - dirbuffer;
            next_id++;

            ((WORD*)&H0)[c] = 0;
            ((WORD*)&F0)[c] = 2 * gap_open_penalty + 2 * gap_extend_penalty;


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
        dprofile_shuffle16(dprofile, score_matrix, dseq);
      else
        dprofile_fill16(dprofile, score_matrix, dseq);

      MQ = v_and(M, Q);
      MR = v_and(M, R);
      MQ0 = MQ;

      domasked16(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
    }

    F0 = v_add_u(F0, R);
    F0 = v_add_u(F0, R);
    F0 = v_add_u(F0, R);
    H0 = v_sub_u(F0, Q);
    F0 = v_add_u(F0, R);


    dir += 4*longestdbsequence;

    if (dir >= dirbuffer + dirbuffersize)
      dir -= dirbuffersize;
  }
}
