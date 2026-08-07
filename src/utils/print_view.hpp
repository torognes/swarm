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

#ifndef SWARM_UTILS_PRINT_VIEW_H
#define SWARM_UTILS_PRINT_VIEW_H


#include "decimal_digits.hpp"  // decimal::Buffer, decimal::to_decimal
#include "view.hpp"  // View<char>
#include <cassert>
#include <cstddef>  // std::size_t
#include <cstdio>  // std::FILE, std::fputc, std::fwrite


// Emit the bytes of a View to a stream.
//
// This is what a "%.*s" conversion was used for before: pass a counted run
// of characters that is not NUL-terminated (a slice of a header inside a
// shared buffer, a decoded sequence). std::fwrite takes both its element
// size and its count as an std::size_t, and the count is what View::size()
// already returns, so nothing is narrowed to the int that "%.*s" requires
// -- and there is no format string to parse at run time.
//
// Not exactly interchangeable with "%.*s": that conversion stops at an
// embedded NUL as well as at the precision, while fwrite always emits the
// full byte count. No header and no sequence reaching swarm's printers
// contains a NUL, because the reader truncates the line at the first one
// (Line_buffer::read_next), so the two agree on every input -- by
// construction of the input rather than by rejection of it. If swarm ever
// rejects NUL instead of truncating, this is the note to re-read.
//
// It lives here rather than in view.hpp so that View stays a pure data type
// and <cstdio> stays out of a header that nearly every translation unit
// reaches through seqinfo.hpp. Same reasoning as view_stream.hpp and
// <ostream>.

inline auto fprint(std::FILE * const output_handle, View<char> const text) -> void
{
  // An empty view may carry a null pointer, and passing one to fwrite is
  // undefined even with a zero count. Reachable: a header whose abundance
  // annotation starts at offset 0 prints an empty leading slice.
  if (text.empty()) { return; }
  // fwrite's two size arguments are an element size and a count, in that
  // order, and transposing them compiles. Naming the element size says
  // which is which; sizeof(char) is 1 by definition, so this is about the
  // reader, not about a platform where it could differ.
  static constexpr std::size_t element_size = sizeof(char);
  static_cast<void>(std::fwrite(text.data(), element_size, text.size(), output_handle));
}


// Emit one character: what std::fputc was used for.
//
// The only thing this adds over the call it wraps is that the discarded
// return value is dealt with once, here, instead of at every call site --
// swarm had 45 'static_cast<void>(std::fputc(...))' spellings. Verified to
// compile to the identical 'jmp fputc' at -O3.
inline auto fprint(std::FILE * const output_handle, char const character) -> void
{
  static_cast<void>(std::fputc(character, output_handle));
}


// Emit a string literal: what std::fputs was used for.
//
// std::fwrite with the array's own bound, not std::fputs, so the length is a
// compile-time constant by construction rather than by optimisation. GCC does
// fold fputs of a literal into exactly this fwrite -- and does so even under
// _FORTIFY_SOURCE, and even through this wrapper, all three checked in the
// disassembly -- but that is a property of one compiler, and spelling it out
// costs nothing. (It is fprintf that fortify stops GCC from folding, which is
// why the fprintf calls elsewhere are worth replacing and these are not.)
//
// The parameter is a reference to an array so that Size arrives with it; a
// char const * would have to be walked at run time. The C array is therefore
// deliberate, as in fatal.hpp's explicit_decay.
//
// Intended for literals, and the contract is the array's bound, not a
// terminator: passing a partially-filled 'char buf[64]' would emit all 63
// bytes, not the string inside it. swarm declares no char arrays, and the
// assert catches the unterminated case; a filled-then-truncated buffer is
// the gap the assert cannot close, hence this note.
// NOLINTNEXTLINE(cppcoreguidelines-avoid-c-arrays,hicpp-avoid-c-arrays,modernize-avoid-c-arrays)
template <std::size_t Size>
auto fprint(std::FILE * const output_handle, char const (&literal)[Size]) -> void
{
  static_assert(Size > 0, "a string literal always carries its terminator");
  assert(literal[Size - 1] == '\0');
  // Size counts the terminating NUL, which is not part of the output. An
  // empty literal leaves count 0, which fwrite accepts: unlike the View
  // overload above, the pointer here cannot be null.
  static constexpr std::size_t element_size = sizeof(char);
  static_cast<void>(std::fwrite(literal, element_size, Size - 1, output_handle));
}


// Emit one integer, in decimal, to a stream: what an "%u" or a "%" PRIu64
// conversion was used for. The digits come from decimal_digits.hpp, so this
// is the same fwrite as above with a locally-produced view.
//
// Deliberately not a batching writer. A tab-separated run of numbers could
// be accumulated in a fixed buffer and emitted with a single fwrite, which
// measures faster in isolation -- but std::FILE * already buffers, so the
// only thing such a buffer saves is stdio call count, i.e. the stream lock.
// Measured over 2.19 M links of a "-d 1 -i" run the whole difference is
// 0.13 % of the runtime, against a second buffer in front of stdio's, its
// zero-initialisation on every construction, and a lifetime to manage at
// every call site. Not worth it; the note is here so the idea is not
// re-proposed. (putc_unlocked would beat both, and is a POSIX-ism that
// MinGW spells differently, so it is out.)
template <typename Integer>
auto fprint_integer(std::FILE * const output_handle, Integer const value) -> void
{
  decimal::Buffer buffer {};
  fprint(output_handle, decimal::to_decimal(buffer, value));
}


// tests:
//
// Four overloads share the name 'fprint', so which one a call picks is the
// part worth pinning down. Checked by capturing the bytes written to a
// std::tmpfile():
//
// fprint(stream, '\t')             -> "\t"        the char overload
// fprint(stream, "\t*\n")          -> "\t*\n"     the literal overload
// fprint(stream, "")               -> ""          count 0, pointer not null
// fprint(stream, View<char>{d, 3}) -> "abc"       the View overload
// fprint(stream, View<char>{})     -> ""          returns before fwrite
// fprint_integer(stream, 4294967295U) -> "4294967295"
//
// No call is ambiguous: a char and a char-array reference are each an exact
// match for their own overload, and View's two-argument constructor is
// explicit, so a literal cannot reach the View overload instead. The
// non-template View overload also outranks the array template for a View
// argument. Verified at -O3 that the char and literal wrappers emit the
// byte-identical instruction stream to the fputc/fputs calls they replace.

#endif // SWARM_UTILS_PRINT_VIEW_H
