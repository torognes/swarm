/*
    SWARM

    Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe

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


#ifndef SWARM_UTILS_BLOOM_H
#define SWARM_UTILS_BLOOM_H

#include <algorithm>  // std::max
#include <atomic>  // std::atomic, std::memory_order_relaxed
#include <cassert>
#include <cstdint>  // uint64_t
#include <limits>  // std::numeric_limits
#include <type_traits>  // std::integral_constant, std::true_type, std::false_type
#include <vector>


namespace bloom_detail {

  // bitmap is stored as an array of 64-bit words; this is the size of
  // one such word (in bytes), used both to lower-bound the requested
  // size and to convert bytes -> words via a right shift.
  constexpr uint64_t bytes_per_word {8};
  static_assert(bytes_per_word == sizeof(uint64_t),
                "bytes_per_word must match sizeof(uint64_t)");

  // The shift that performs that bytes -> words conversion. Named, and
  // checked against bytes_per_word, because the constructor used to spell
  // it as a bare 3: the comment above claimed bytes_per_word was what
  // converted, while nothing tied the two together.
  constexpr unsigned int bytes_per_word_shift {3};
  static_assert((uint64_t{1} << bytes_per_word_shift) == bytes_per_word,
                "bytes_per_word_shift must halve as many times as bytes_per_word divides");

  // Give each pattern n_hash_functions distinct bits, drawn at random.
  //
  // Not a member and not a template: the pattern *shift* decides how many
  // patterns there are, but nothing here depends on it, so this loop is
  // compiled once for the whole program instead of once per instantiation
  // -- and <random> stays out of a header. Defined in bloom.cpp.
  //
  // Not marked noexcept: rand_64.operator() (std::mt19937_64) is not
  // formally noexcept in the standard, even if it does not throw in
  // practice. Called only from the constructor, which is itself non-
  // noexcept (vector resize/construction can throw bad_alloc), so the
  // distinction is academic.
  auto generate_patterns(std::vector<uint64_t> & patterns,
                         uint64_t pattern_k) -> void;

}  // namespace bloom_detail


// How a filter's bitmap length relates to a power of two, which decides how
// a hash is mapped onto it (see bitmap_index). This is a property of where
// the length comes from, and every construction site knows it:
// compute_hashtable_size() returns an exact power of two, the --bloom-bits /
// --ceiling budget returns whatever it returns.
enum struct Bitmap_size : std::uint8_t { power_of_two, arbitrary };


// Blocked Bloom filter with precomputed bit patterns
// (Putze, Sanders, Singler 2009 -- see bloom.cpp for the reference).
// Bit semantics are inverted from a textbook Bloom filter: a freshly
// constructed filter has all bits set to 1, set() clears the pattern's
// bits, and get() returns true (possibly-present) when all of the
// pattern's bits in the addressed word are zero.
//
// pattern_shift is a template parameter rather than a constructor one:
// every filter swarm builds names it with a constant (amplicon_pattern_shift
// and fastidious_pattern_shift, both in algod1_internal.hpp), so the
// runtime parameter only ever carried a compile-time fact. Lifting it into
// the type turns the shift in bitmap_index() and the mask in bit_pattern()
// from loaded members into immediates, and removes four of the five scalar
// members: the shift itself, the count and mask derived from it, and
// pattern_k, which only the pattern generation above ever read. Both are on
// the probe path, which runs twice per Bloom lookup in the fastidious
// passes.
//
// It also gives the two filters distinct types. Six functions in
// algod1_fastidious.cpp take bloom_a and bloom_f as adjacent parameters,
// which were the same type and so could be passed in either order; now
// they cannot.
//
// n_hash_functions stays a constructor parameter: the fastidious filter
// derives its k from --bloom-bits (see compute_bloom_geometry), and k is
// read only by the cold pattern generation, so there would be nothing to
// fold.
template <unsigned int pattern_shift, Bitmap_size size_kind>
class BloomFilter {
public:

  // bitmap_bytes is the requested bitmap size in bytes; it is rounded
  // up to at least one 64-bit word so bitmap_index() can compute a
  // valid position.
  //
  // Non-noexcept: the two vector constructions can throw std::bad_alloc.
  //
  // The words start all-ones through an explicit store loop rather than
  // the old (count, value) construction: std::atomic is not copyable, so
  // the vector can only be built default-initialized. Construction is
  // single-threaded, but the loop stays on the atomic interface -- mixing
  // atomic and plain access to the same object is undefined behaviour.
  BloomFilter(uint64_t const bitmap_bytes,
              unsigned int const n_hash_functions)
    : size{std::max(bitmap_bytes, bloom_detail::bytes_per_word)
            >> bloom_detail::bytes_per_word_shift}
    , bitmap(size)
    , patterns(pattern_count) {
    // Checked here, once, rather than at every probe: this is the promise
    // the caller makes by choosing a size_kind, and the constructor is
    // where the value it is about arrives. A power-of-two filter whose
    // length is not one would silently address a fraction of its bitmap --
    // see bitmap_index.
    assert(size != 0);
    assert(size_kind == Bitmap_size::arbitrary or (size & (size - 1)) == 0);
    for (auto & word : bitmap) {
      word.store(std::numeric_limits<uint64_t>::max(),
                 std::memory_order_relaxed);
    }
    bloom_detail::generate_patterns(patterns, n_hash_functions);
  }

  // Mark hash as a member of the set.
  //
  // An atomic fetch_and rather than a plain '&=': the fastidious light
  // pass calls set() from every worker thread with no lock held
  // (mark_light_thread in algod1_fastidious.cpp), and two plain
  // read-modify-writes landing on the same word could each lose the
  // other's bits -- a lost pattern turns get() into a false *negative*,
  // which is the one error a Bloom filter must never make (a missed
  // graft, and a run whose output depends on thread timing). Relaxed
  // ordering suffices: concurrent set() calls need atomicity only, and
  // every reader runs after the writers' ThreadRunner has joined, which
  // already orders the passes.
  auto set(uint64_t const hash) noexcept -> void {
    bitmap[bitmap_index(hash)].fetch_and(compl bit_pattern(hash),
                                         std::memory_order_relaxed);
  }

  // Test whether hash may be a member of the set. Returns true on
  // possible-membership, false on definite-non-membership.
  auto get(uint64_t const hash) const noexcept -> bool {
    return (bitmap[bitmap_index(hash)].load(std::memory_order_relaxed)
            & bit_pattern(hash)) == 0U;
  }

private:

  // Derived from the shift, so constants rather than members.
  static constexpr uint64_t pattern_count {uint64_t{1} << pattern_shift};
  static constexpr uint64_t pattern_mask {pattern_count - 1};

  // The largest bitmap the multiply-shift map below can address: its
  // product must fit in a uint64_t, and one factor is already 32 bits
  // wide. 2^32 words is a 32 GiB bitmap.
  static constexpr uint64_t max_mapped_size {uint64_t{1} << 32U};

  // Which word of the bitmap a hash addresses.
  //
  // This used to be '% size' for every filter, and a 64-bit divide is 20 to
  // 40 cycles that cannot be hoisted -- the dividend changes every call --
  // on a path that runs once per microvariant, some 200 million times on a
  // 100k-record 'd = 1 -f' run. Both forms below replace it, and which one
  // applies is decided by the type rather than by a member, so neither
  // filter pays for the other's generality.
  //
  // Note for anyone changing this function: swarm's output cannot tell you
  // whether you got it right. Every candidate the filter passes is verified
  // (check_variant, hash_check_attach), so the mapping only moves the
  // false-positive set; as long as set() and get() agree -- they call this
  // one function, so they do -- the clustering is identical however bad the
  // map is. Masking a non-power-of-two size was measured at 31x slower with
  // byte-identical output. The check that matters here is the clock.
  using Size_is_power_of_two =
    std::integral_constant<bool, size_kind == Bitmap_size::power_of_two>;

  // A power of two: '& (size - 1)' is '% size', exactly. This is the
  // amplicon filter, whose length is compute_hashtable_size()'s return
  // value; the constructor asserts the invariant this relies on.
  auto bitmap_index(uint64_t const hash,
                    std::true_type /* size is a power of two */) const noexcept -> uint64_t {
    auto const position = (hash >> pattern_shift) & (size - 1);
    assert(position < bitmap.size());
    return position;
  }

  // Any length: multiply-shift, which maps the top 32 bits of the hash
  // uniformly onto [0, size) as the high half of a 32x32 product. One
  // multiply, no rounding of the caller's memory budget, and it needs no
  // modulo because the filter never needed one -- it needs a uniform
  // address, and any uniform map is as good as another here. Bucket sizes
  // differ by at most one part in 2^32/size, which would matter to a
  // partitioning scheme and does not matter to a Bloom filter.
  //
  // The bits stay disjoint from bit_pattern's, which reads the low
  // pattern_shift bits, so the two remain independent as they were when
  // this shifted by pattern_shift first.
  //
  // Beyond 32 GiB the product would overflow and the modulo comes back.
  // That path is reachable -- the bitmap is sized from the light clusters'
  // total length, and swarm aims at tera-scale input -- so it is a branch
  // and not an assertion. It is decided by a member that never changes, so
  // it predicts perfectly, and it costs a compare where the divide it
  // guards costs tens of cycles. A 128-bit product would remove even that,
  // but __int128 is not ISO C++ and -Wpedantic is one of this project's
  // warnings.
  auto bitmap_index(uint64_t const hash,
                    std::false_type /* size is arbitrary */) const noexcept -> uint64_t {
    auto const position = (size < max_mapped_size)
      ? (((hash >> 32U) * size) >> 32U)
      : ((hash >> pattern_shift) % size);
    assert(position < bitmap.size());
    return position;
  }

  auto bitmap_index(uint64_t const hash) const noexcept -> uint64_t {
    return bitmap_index(hash, Size_is_power_of_two{});
  }

  auto bit_pattern(uint64_t const hash) const noexcept -> uint64_t {
    auto const position = hash & pattern_mask;
    assert(position < patterns.size());
    return patterns[position];
  }

  uint64_t size {0};            // bitmap length, in 64-bit words
  // atomic words: see set() -- same size and layout as plain uint64_t,
  // and lock-free on every supported target (x86_64, aarch64, ppc64le,
  // x86_64-mingw)
  std::vector<std::atomic<uint64_t>> bitmap;
  std::vector<uint64_t> patterns;
};

#endif  // SWARM_UTILS_BLOOM_H
