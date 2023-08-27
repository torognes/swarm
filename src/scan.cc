/*
    SWARM

    Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe

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
#include "utils/nt_codec.h"
#include <cassert>  // assert()
#include <climits>
#include <cstdint>  // int64_t, uint64_t
#include <pthread.h>


static pthread_mutex_t scan_mutex;

static ThreadRunner * search_threads = nullptr;

struct Search_data
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

static struct Search_data * sd;
static uint64_t master_next;
static uint64_t master_length;
static uint64_t remainingchunks;
static uint64_t * master_targets;
static uint64_t * master_scores;
static uint64_t * master_diffs;
static uint64_t * master_alignlengths;
static int master_bits;
static uint64_t dirbuffersize;

queryinfo_t query;
uint64_t longestdbsequence;


void search_alloc(struct Search_data * search_data)
{
  static constexpr unsigned int one_kilobyte {1024};
  static constexpr unsigned int nt_per_uint64 {32};

  dirbuffersize = longestdbsequence * ((longestdbsequence + 3) / 4) * 4;
  search_data->qtable = new BYTE*[longestdbsequence];
  search_data->qtable_w = new WORD*[longestdbsequence];
  search_data->dprofile = new BYTE[2 * one_kilobyte];  // 4 * 16 * 32
  search_data->dprofile_w = new WORD[1 * one_kilobyte];  // 4 * 2 * 8 * 32
  search_data->hearray = new BYTE[longestdbsequence * nt_per_uint64] { };
  search_data->dir_array = new uint64_t[dirbuffersize] { };
}

void search_free(struct Search_data * search_data)
{
  delete [] search_data->qtable;
  search_data->qtable = nullptr;
  delete [] search_data->qtable_w;
  search_data->qtable_w = nullptr;
  delete [] search_data->dprofile;
  search_data->dprofile = nullptr;
  delete [] search_data->dprofile_w;
  search_data->dprofile_w = nullptr;
  delete [] search_data->hearray;
  search_data->hearray = nullptr;
  delete [] search_data->dir_array;
  search_data->dir_array = nullptr;
}

void search_init(struct Search_data * search_data)
{
  static constexpr int byte_multiplier {64};
  static constexpr int word_multiplier {32};

  for(auto i = 0U; i < query.len; i++)
  {
    const int nt_value {nt_extract(query.seq, i) + 1};   // 1,  2,   3, or   4
    const int byte_offset {byte_multiplier * nt_value};  // 1, 64, 128, or 192
    const int word_offset {word_multiplier * nt_value};  // 1, 32,  64, or 128

    search_data->qtable[i]   = search_data->dprofile   + byte_offset;
    search_data->qtable_w[i] = search_data->dprofile_w + word_offset;
  }
}

void search_chunk(struct Search_data * search_data, const int64_t bits)
{
  static constexpr unsigned int bit_mode_16 {16};

  assert(search_data->target_count != 0);
  assert((bits == bit_mode_16) || (bits == bit_mode_16 / 2));

 if (bits == bit_mode_16)
   {
    search16(search_data->qtable_w,
             static_cast<WORD>(penalty_gapopen),
             static_cast<WORD>(penalty_gapextend),
             static_cast<WORD*>(score_matrix_16),
             search_data->dprofile_w,
             reinterpret_cast<WORD*>(search_data->hearray),
             search_data->target_count,
             master_targets + search_data->target_index,
             master_scores + search_data->target_index,
             master_diffs + search_data->target_index,
             master_alignlengths + search_data->target_index,
             static_cast<uint64_t>(query.len),
             dirbuffersize,
             search_data->dir_array,
             longestdbsequence);
   }
 else {
    search8(search_data->qtable,
            static_cast<BYTE>(penalty_gapopen),
            static_cast<BYTE>(penalty_gapextend),
            static_cast<BYTE*>(score_matrix_8),
            search_data->dprofile,
            search_data->hearray,
            search_data->target_count,
            master_targets + search_data->target_index,
            master_scores + search_data->target_index,
            master_diffs + search_data->target_index,
            master_alignlengths + search_data->target_index,
            static_cast<uint64_t>(query.len),
            dirbuffersize,
            search_data->dir_array,
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
      const uint64_t chunksize =
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


void search_worker_core(const int64_t t)
{
  search_init(sd + t);
  while(search_getwork(& sd[t].target_count, & sd[t].target_index)) {
    search_chunk(sd + t, master_bits);
  }
}


auto adjust_thread_number(const int n_bits,
                          const uint64_t remaining_sequences,
                          uint64_t n_threads) -> uint64_t {
  static constexpr unsigned int channels_8 {8};
  static constexpr unsigned int channels_16 {16};
  static constexpr unsigned int bit_mode_16 {16};  // refactoring: should be an enum class
  const auto channels = (n_bits == bit_mode_16) ? channels_8 : channels_16;

  assert(remaining_sequences != 0);
  assert(n_threads != 0);
  assert((n_threads - 1) <= (ULLONG_MAX / channels_8));
  assert((n_bits == bit_mode_16) || (n_bits == bit_mode_16 / 2));

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
               const int bits)
{
  unsigned int query_len = 0;
  query.qno = query_no;
  db_getsequenceandlength(query_no, &query.seq, &query_len);
  query.len = query_len;

  master_next = 0;
  master_length = listlength;
  master_targets = targets;
  master_scores = scores;
  master_diffs = diffs;
  master_alignlengths = alignlengths;
  master_bits = bits;

  const auto thr =
    adjust_thread_number(bits,
                         master_length,
                         static_cast<uint64_t>(opt_threads));

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

  sd = new struct Search_data[static_cast<uint64_t>(opt_threads)];

  for(auto t = 0LL; t < opt_threads; t++) {
    search_alloc(sd + t);
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
    search_free(sd + t);
  }
  delete [] sd;
  sd = nullptr;
}
