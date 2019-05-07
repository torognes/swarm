
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

#ifdef __aarch64__

typedef int16x8_t VECTORTYPE;

const uint16x8_t neon_mask =
  {0x0003, 0x000c, 0x0030, 0x00c0, 0x0300, 0x0c00, 0x3000, 0xc000};

const VECTORTYPE T0_init = { -1, 0, 0, 0, 0, 0, 0, 0 };

#define v_init(a,b,c,d,e,f,g,h) (const VECTORTYPE){a,b,c,d,e,f,g,h}
#define v_load(a) vld1q_s16((const int16_t *)(a))
#define v_store(a, b) vst1q_s16((int16_t *)(a), (b))
#define v_merge_lo_16(a, b) vzip1q_s16((a),(b))
#define v_merge_hi_16(a, b) vzip2q_s16((a),(b))
#define v_merge_lo_32(a, b) vreinterpretq_s16_s32(vzip1q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)))
#define v_merge_hi_32(a, b) vreinterpretq_s16_s32(vzip2q_s32(vreinterpretq_s32_s16(a), vreinterpretq_s32_s16(b)))
#define v_merge_lo_64(a, b) vreinterpretq_s16_s64(vcombine_s64(vget_low_s64(vreinterpretq_s64_s16(a)), vget_low_s64(vreinterpretq_s64_s16(b))))
#define v_merge_hi_64(a, b) vreinterpretq_s16_s64(vcombine_s64(vget_high_s64(vreinterpretq_s64_s16(a)), vget_high_s64(vreinterpretq_s64_s16(b))))
#define v_min(a, b) vminq_s16((a), (b))
#define v_min_u(a, b) vminq_u16((a), (b))
#define v_add(a, b) vqaddq_u16((a), (b))
#define v_sub(a, b) vqsubq_u16((a), (b))
#define v_dup(a) vdupq_n_s16(a)
#define v_zero v_dup(0)
#define v_and(a, b) vandq_u16((a), (b))
#define v_xor(a, b) veorq_u16((a), (b))
#define v_shift_left(a) vextq_u16((v_zero), (a), 7)
#define v_mask_gt(a, b) vaddvq_u16(vandq_u16((vcgtq_s16((a), (b))), neon_mask))
#define v_mask_eq(a, b) vaddvq_u16(vandq_u16((vceqq_s16((a), (b))), neon_mask))

#elif __x86_64__

typedef __m128i VECTORTYPE;

const VECTORTYPE T0_init = _mm_set_epi16(0, 0, 0, 0, 0, 0, 0, -1);

#define v_load(a) _mm_load_si128((VECTORTYPE *)(a))
#define v_store(a, b) _mm_store_si128((VECTORTYPE *)(a), (b))
#define v_merge_lo_16(a, b) _mm_unpacklo_epi16((a),(b))
#define v_merge_hi_16(a, b) _mm_unpackhi_epi16((a),(b))
#define v_merge_lo_32(a, b) _mm_unpacklo_epi32((a),(b))
#define v_merge_hi_32(a, b) _mm_unpackhi_epi32((a),(b))
#define v_merge_lo_64(a, b) _mm_unpacklo_epi64((a),(b))
#define v_merge_hi_64(a, b) _mm_unpackhi_epi64((a),(b))
#define v_min(a, b) _mm_min_epi16((a), (b))
#define v_min_u(a, b) _mm_min_epu16((a), (b))
#define v_add(a, b) _mm_adds_epu16((a), (b))
#define v_sub(a, b) _mm_subs_epu16((a), (b))
#define v_dup(a) _mm_set1_epi16(a)
#define v_zero v_dup(0)
#define v_and(a, b) _mm_and_si128((a), (b))
#define v_xor(a, b) _mm_xor_si128((a), (b))
#define v_shift_left(a) _mm_slli_si128((a), 2)
#define v_mask_gt(a, b) _mm_movemask_epi8(_mm_cmpgt_epi16((a), (b)))
#define v_mask_eq(a, b) _mm_movemask_epi8(_mm_cmpeq_epi16((a), (b)))

#elif __PPC__

#error Architecture ppcle64 not implemented yet

#else

#error Unknown Architecture

#endif

void dprofile_dump16(WORD * dprofile)
{
  printf("\ndprofile:\n");
  for(int i=0; i<32; i++)
    {
      printf("%c: ", sym_nt[i]);
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

inline void onestep_16(VECTORTYPE & H,
                       VECTORTYPE & N,
                       VECTORTYPE & F,
                       VECTORTYPE V,
                       unsigned short * DIR,
                       VECTORTYPE & E,
                       VECTORTYPE QR,
                       VECTORTYPE R)
{
  H = v_add(H, V);
  *(DIR+0) = v_mask_gt(H, F);
  H = v_min(H, F);
  H = v_min(H, E);
  *(DIR+1) = v_mask_eq(H, E);
  N = H;
  H = v_add(H, QR);
  F = v_add(F, R);
  E = v_add(E, R);
  *(DIR+2) = v_mask_gt(H, F);
  *(DIR+3) = v_mask_gt(H, E);
  F = v_min(H, F);
  E = v_min(H, E);
}

void align_cells_regular_16(VECTORTYPE * Sm,
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
  h5 = v_zero;
  h6 = v_zero;
  h7 = v_zero;
  h8 = v_zero;

  for(long i = 0; i < ql; i++)
    {
      x = qp[i + 0];
      h4 = hep[2*i + 0];
      E  = hep[2*i + 1];
      onestep_16(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R);
      onestep_16(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R);
      onestep_16(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R);
      onestep_16(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R);
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

void align_cells_masked_16(VECTORTYPE * Sm,
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
  h5 = v_zero;
  h6 = v_zero;
  h7 = v_zero;
  h8 = v_zero;

  for(long i = 0; i < ql; i++)
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

      onestep_16(h0, h5, f0, x[0], dir + 16*i +  0, E, Q, R);
      onestep_16(h1, h6, f1, x[1], dir + 16*i +  4, E, Q, R);
      onestep_16(h2, h7, f2, x[2], dir + 16*i +  8, E, Q, R);
      onestep_16(h3, h8, f3, x[3], dir + 16*i + 12, E, Q, R);
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

unsigned long backtrack_16(char * qseq,
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

      unsigned long d = dirbuffer[(offset + longestdbsequence*4*(j/4)
                                   + 4*i + (j&3)) % dirbuffersize];

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

  T0 = T0_init;
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

          align_cells_regular_16(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0);
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
                      long score = ((WORD*)S)[z*CHANNELS+c];
                      scores[cand_id] = score;

                      unsigned long diff;

                      if (score < 65535)
                        {
                          long offset = d_offset[c];
                          diff = backtrack_16(query.seq, dbseq, qlen, dbseqlen,
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

          align_cells_masked_16(S, hep, qp, &Q, &R, qlen, &F0, dir, &H0, &M, &MQ, &MR, &MQ0);
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
