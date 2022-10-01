/*
    SWARM

    Copyright (C) 2012-2022 Torbjorn Rognes and Frederic Mahe

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
#include "db.h"
#include "search8.h"
#include "search16.h"
#include "threads.h"


static pthread_mutex_t scan_mutex;

static ThreadRunner * search_threads = nullptr;

struct search_data
{
  BYTE ** qtable;
  WORD ** qtable_w;

  BYTE * dprofile;
  WORD * dprofile_w;

  BYTE * hearray;

  uint64_t * dir_array;

  uint64_t target_count;
  uint64_t target_index;
};

static struct search_data * sd;
static uint64_t master_next;
static uint64_t master_length;
static uint64_t remainingchunks;
static uint64_t * master_targets;
static uint64_t * master_scores;
static uint64_t * master_diffs;
static uint64_t * master_alignlengths;
static int master_bits;
static uint64_t dirbufferbytes;

queryinfo_t query;
uint64_t longestdbsequence;


void search_alloc(struct search_data * sdp)
{
  constexpr unsigned int one_kilobyte {1024};
  constexpr unsigned int nt_per_uint64 {32};
  constexpr unsigned int bytes_per_uint64 {8};

  dirbufferbytes =
    bytes_per_uint64 * longestdbsequence * ((longestdbsequence + 3) / 4) * 4;
  sdp->qtable = new BYTE*[longestdbsequence];
  sdp->qtable_w = new WORD*[longestdbsequence];
  sdp->dprofile = new BYTE[2 * one_kilobyte];  // 4 * 16 * 32
  sdp->dprofile_w = new WORD[2 * one_kilobyte];  // 4 * 2 * 8 * 32
  sdp->hearray = new BYTE[longestdbsequence * nt_per_uint64] { };
  sdp->dir_array = new uint64_t[dirbufferbytes] { };
}

void search_free(struct search_data * sdp)
{
  delete [] sdp->qtable;
  sdp->qtable = nullptr;
  delete [] sdp->qtable_w;
  sdp->qtable_w = nullptr;
  delete [] sdp->dprofile;
  sdp->dprofile = nullptr;
  delete [] sdp->dprofile_w;
  sdp->dprofile_w = nullptr;
  delete [] sdp->hearray;
  sdp->hearray = nullptr;
  delete [] sdp->dir_array;
  sdp->dir_array = nullptr;
}

void search_init(struct search_data * sdp)
{
  constexpr int byte_multiplier {64};
  constexpr int word_multiplier {32};

  for(auto i = 0U; i < query.len; i++)
  {
    const int nt_value {nt_extract(query.seq, i) + 1};   // 1,  2,   3, or   4
    const int byte_offset {byte_multiplier * nt_value};  // 1, 64, 128, or 192
    const int word_offset {word_multiplier * nt_value};  // 1, 32,  64, or 128

    sdp->qtable[i]   = sdp->dprofile   + byte_offset;
    sdp->qtable_w[i] = sdp->dprofile_w + word_offset;
  }
}

void search_chunk(struct search_data * sdp, const int64_t bits)
{
  constexpr unsigned int bit_mode_16 {16};
  if (sdp->target_count == 0) {
    return;
  }

 if (bits == bit_mode_16)
   {
    search16(sdp->qtable_w,
             static_cast<WORD>(penalty_gapopen),
             static_cast<WORD>(penalty_gapextend),
             static_cast<WORD*>(score_matrix_16),
             sdp->dprofile_w,
             reinterpret_cast<WORD*>(sdp->hearray),
             sdp->target_count,
             master_targets + sdp->target_index,
             master_scores + sdp->target_index,
             master_diffs + sdp->target_index,
             master_alignlengths + sdp->target_index,
             static_cast<uint64_t>(query.len),
             dirbufferbytes / sizeof(uint64_t),
             sdp->dir_array,
             longestdbsequence);
   }
 else {
    search8(sdp->qtable,
            static_cast<BYTE>(penalty_gapopen),
            static_cast<BYTE>(penalty_gapextend),
            static_cast<BYTE*>(score_matrix_8),
            sdp->dprofile,
            sdp->hearray,
            sdp->target_count,
            master_targets + sdp->target_index,
            master_scores + sdp->target_index,
            master_diffs + sdp->target_index,
            master_alignlengths + sdp->target_index,
            static_cast<uint64_t>(query.len),
            dirbufferbytes / sizeof(uint64_t),
            sdp->dir_array,
            longestdbsequence);
 }
}

auto search_getwork(uint64_t * countref, uint64_t * firstref) -> bool
{
  // * countref = how many sequences to search
  // * firstref = index into master_targets/scores/diffs where thread should start

  bool status {false};

  pthread_mutex_lock(&scan_mutex);

  if (master_next < master_length)
    {
      uint64_t chunksize =
        ((master_length - master_next + remainingchunks - 1) / remainingchunks);

      * countref = chunksize;
      * firstref = master_next;

      master_next += chunksize;
      remainingchunks--;
      status = true;
    }

  pthread_mutex_unlock(&scan_mutex);

  return status;
}

void search_worker_core(int64_t t)
{
  search_init(sd + t);
  while(search_getwork(& sd[t].target_count, & sd[t].target_index)) {
    search_chunk(sd + t, master_bits);
  }
}


auto adjust_thread_number(const int n_bits,
                          const uint64_t remaining_sequences,
                          uint64_t n_threads) -> uint64_t {
  constexpr unsigned int channels_8 {8};
  constexpr unsigned int channels_16 {16};
  constexpr unsigned int bit_mode_16 {16};
  auto channels {channels_16};

  assert(remaining_sequences > 0);
  assert(n_threads > 0);
  assert((n_bits == bit_mode_16) || (n_bits == bit_mode_16 / 2));

  if (n_bits == bit_mode_16) {
    channels = channels_8;
  }

  while (remaining_sequences <= (n_threads - 1) * channels) {
    --n_threads;
  }

  return n_threads;
}

// arguments: bits, master_length, thr
// static_assert(adjust_thread_number( 8, 32, 10) == 2);
// static_assert(adjust_thread_number( 8, 32,  3) == 2);
// static_assert(adjust_thread_number( 8, 31,  2) == 2);
// static_assert(adjust_thread_number( 8, 17,  2) == 2);
// static_assert(adjust_thread_number( 8, 16,  2) == 1);
// static_assert(adjust_thread_number( 8,  1,  2) == 1);
// static_assert(adjust_thread_number( 8, 32,  1) == 1);
// static_assert(adjust_thread_number(16, 17, 10) == 3);
// static_assert(adjust_thread_number(16, 17,  3) == 3);
// static_assert(adjust_thread_number(16, 16,  3) == 2);
// static_assert(adjust_thread_number(16, 15,  2) == 2);
// static_assert(adjust_thread_number(16,  1,  3) == 1);
// static_assert(adjust_thread_number(16, 17,  1) == 1);


void search_do(uint64_t query_no,
               uint64_t listlength,
               uint64_t * targets,
               uint64_t * scores,
               uint64_t * diffs,
               uint64_t * alignlengths,
               int bits)
{
  // constexpr unsigned int bit_mode_16 {16};
  query.qno = query_no;
  unsigned int query_len = 0;
  db_getsequenceandlength(query_no, &query.seq, &query_len);
  query.len = query_len;

  master_next = 0;
  master_length = listlength;
  master_targets = targets;
  master_scores = scores;
  master_diffs = diffs;
  master_alignlengths = alignlengths;
  master_bits = bits;

  auto thr = static_cast<uint64_t>(opt_threads);

  thr = adjust_thread_number(bits, master_length, thr);

  remainingchunks = thr;

  if (thr == 1) {
    search_worker_core(0);
  }
  else {
    search_threads->run();
  }
}

void search_begin()
{
  longestdbsequence = db_getlongestsequence();

  sd = new struct search_data[static_cast<uint64_t>(opt_threads)];

  for(auto t = 0LL; t < opt_threads; t++) {
    search_alloc(sd+t);
  }

  pthread_mutex_init(& scan_mutex, nullptr);

  /* start threads */

  search_threads
    = new ThreadRunner(static_cast<int>(opt_threads), search_worker_core);
}

void search_end()
{
  /* finish and clean up worker threads */

  delete search_threads;
  search_threads = nullptr;

  pthread_mutex_destroy(& scan_mutex);

  for(auto t = 0LL; t < opt_threads; t++) {
    search_free(sd+t);
  }
  delete [] sd;
  sd = nullptr;
}
