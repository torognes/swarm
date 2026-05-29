# SWARM
#
# Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
# Department of Informatics, University of Oslo,
# PO Box 1080 Blindern, NO-0316 Oslo, Norway

# Makefile for SWARM

COMMON = -g -fno-exceptions -pthread
LINKOPT =
LIBS =
WARNINGS = -Wall -Wextra -Wpedantic
# Generate header-dependency files alongside each object so that
# editing a header forces a rebuild of the right translation units.
DEPFLAGS = -MMD -MP

# Extra warnings used in DEBUG builds. Split per compiler family
# because GCC and clang each understand a slightly different subset.
DEBUG_WARNINGS_COMMON = -Wcast-align -Wcast-qual -Wconversion -Wdate-time \
	-Wdouble-promotion -Wfloat-equal -Wformat=2 -Wnon-virtual-dtor \
	-Wnull-dereference -Wold-style-cast -Woverloaded-virtual -Wshadow \
	-Wsign-conversion -Wuninitialized -Wunused -Wunused-macros -Wvla
DEBUG_WARNINGS_GCC = -Wduplicated-branches -Wduplicated-cond -Wformat-overflow \
	-Wlogical-op -Wunsafe-loop-optimizations -Wuseless-cast
DEBUG_WARNINGS_CLANG = -Wextra-semi -Wcomma -Wassign-enum -Wover-aligned
DEBUG_SANITIZER = -fsanitize=undefined,address -fno-omit-frame-pointer


# Compiler identity and major version. Use the preprocessor to detect
# clang via the __clang__ macro: that works regardless of the binary
# name (g++, clang++, ccache wrappers, cross-compilers...). Use $(CXX)
# -dumpversion to query the actual compiler being used, not the host
# /usr/bin/gcc.
IS_CLANG := $(shell $(CXX) -x c++ -E -dM - < /dev/null 2>/dev/null | grep -c '__clang__')
# Some toolchains report a suffixed major version (e.g. mingw's gcc
# emits '13-win32'); keep only the leading digit run so the numeric
# comparison further down ('-ge 12') stays well-formed.
CXX_VERSION_MAJOR := $(shell $(CXX) -dumpversion 2>/dev/null | sed -E 's/[^0-9].*//')
CXX_VERSION_MAJOR := $(or $(CXX_VERSION_MAJOR),0)

# Pick the right extra-warning set for the detected compiler.
ifneq ($(IS_CLANG), 0)
	DEBUG_WARNINGS = $(DEBUG_WARNINGS_COMMON) $(DEBUG_WARNINGS_CLANG)
else
	DEBUG_WARNINGS = $(DEBUG_WARNINGS_COMMON) $(DEBUG_WARNINGS_GCC)
endif


# Default to RELEASE when no other build flavour was requested. Doing
# this here lets us write the RELEASE flags exactly once below.
ifeq ($(or $(RELEASE),$(PROFILE),$(DEBUG),$(COVERAGE)),)
	RELEASE := 1
endif

# Run "make RELEASE=1" to compile for release (also the default)
ifdef RELEASE
	COMMON += -O3 -DNDEBUG
	ifneq ($(IS_CLANG), 0)
		COMMON += -flto=full
	else
		ifeq ($(shell [ $(CXX_VERSION_MAJOR) -ge 12 ] && echo yes), yes)
			COMMON += -flto=auto
		else
			COMMON += -flto
		endif
	endif

# Run "make PROFILE=1" to compile for profiling
else ifdef PROFILE
	COMMON += -pg -O1

# Run "make DEBUG=1" to compile for debugging (GCC 10 or more recent)
else ifdef DEBUG
	COMMON += -O0 -DDEBUG -ggdb3 -D_GLIBCXX_DEBUG $(DEBUG_SANITIZER) $(DEBUG_WARNINGS)

# Run "make COVERAGE=1" to compile for coverage tests
else ifdef COVERAGE
	COMMON += -fprofile-arcs -ftest-coverage -O0
	LINKOPT += --coverage
	LIBS += -lgcov
endif


# Identify the target machine and OS by asking the compiler. This is
# robust to ccache wrappers, absolute paths, version-suffixed binaries
# (e.g. clang++-19), and clang's --target= cross-compilation.
TARGET_TRIPLE := $(shell $(CXX) -dumpmachine 2>/dev/null)
MACHINE := $(firstword $(subst -, ,$(TARGET_TRIPLE)))
ifeq ($(MACHINE),)
	MACHINE := $(shell uname -m)
endif
# normalise the GCC name "powerpc64le" to the short form used below
ifeq ($(MACHINE), powerpc64le)
	MACHINE := ppc64le
endif


# Machine specific
#
# ARCH_DIR points at the per-ISA folder under src/arch/ that supplies the
# vector intrinsics and any other arch-only translation units. Every .cc
# in that folder gets compiled into the binary (see ARCH_OBJS below).
ifeq ($(MACHINE), x86_64)
	COMMON += -march=x86-64 -mtune=generic -std=c++11
	ARCH_DIR := arch/x86_64
else ifeq ($(MACHINE), aarch64)
	COMMON += -march=armv8-a+simd -mtune=generic \
	          -flax-vector-conversions -std=c++11
	ARCH_DIR := arch/aarch64
else ifeq ($(MACHINE), arm64)
	COMMON += -std=c++11
	ARCH_DIR := arch/aarch64
else ifeq ($(MACHINE), ppc64le)
	COMMON += -mcpu=power8 -maltivec -std=gnu++11
	ARCH_DIR := arch/ppc
endif

# Auto-collect every .cc under the chosen ARCH_DIR. Per-flag overrides
# (e.g. -mssse3 / -msse4.1 / -mpopcnt for the x86_64 dispatch helpers)
# are attached as explicit target rules further down — those win over the
# generic %.o pattern rule.
ARCH_OBJS := $(patsubst %.cc, %.o, $(wildcard $(ARCH_DIR)/*.cc))


# Opt-in portable qgram fallback. Replaces the per-ISA qgram_compare.o
# with arch/portable/qgram_compare.o, a scalar __builtin_popcountll
# loop. Intended as a building block when porting swarm to an
# architecture that has no hand-written SIMD path yet — does NOT
# replace the search8/search16 intrinsics, which still need per-ISA
# implementations.
ifdef SWARM_PORTABLE_FALLBACK
	ARCH_OBJS := $(filter-out $(ARCH_DIR)/qgram_compare.o, $(ARCH_OBJS))
	ARCH_OBJS += arch/portable/qgram_compare.o
endif


# Machine OS specific
#
# OS_DIR points at the per-OS folder under src/os/ that supplies the
# implementations of host-OS-specific helpers (e.g. system_memory.cc).
# The TARGET_TRIPLE substring (linux / mingw / darwin) selects the
# right folder; every .cc in that folder is compiled into the binary
# (see OS_OBJS below). Mirrors the ARCH_DIR pattern used above.
ifneq (,$(findstring linux, $(TARGET_TRIPLE)))
	OS_DIR := os/linux
	BIN = swarm
else ifneq (,$(findstring mingw, $(TARGET_TRIPLE)))
	OS_DIR := os/windows
	LIBS += -lpsapi
	LINKOPT += -static
	BIN = swarm.exe
else ifneq (,$(findstring darwin, $(TARGET_TRIPLE)))
	OS_DIR := os/macos
	BIN = swarm
else ifneq (,$(findstring freebsd, $(TARGET_TRIPLE)))
	OS_DIR := os/freebsd
	BIN = swarm
else
	BIN = swarm
endif

OS_OBJS := $(patsubst %.cc, %.o, $(wildcard $(OS_DIR)/*.cc))


LINKFLAGS = $(COMMON) $(LINKOPT)

CXXFLAGS = $(COMMON) $(WARNINGS)

# Local link target ($(BIN) varies between "swarm" and "swarm.exe")
PROG := $(BIN)
# Final, installed location of the binary, relative to src/
INSTALLED := ../bin/$(BIN)

OBJS = algo.o algod1.o cli.o db.o dereplicate.o \
	qgram.o scan.o search16.o search8.o \
	swarm.o variants.o \
	$(patsubst %.cc, %.o, $(wildcard utils/*.cc)) $(ARCH_OBJS) $(OS_OBJS)

DEPFILES = $(OBJS:.o=.d)

.PHONY: all clean

all : $(INSTALLED)

$(PROG) : $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

# Order-only prerequisite on ../bin: don't relink just because the
# directory's mtime changed.
$(INSTALLED) : $(PROG) | ../bin
	cp -a $(PROG) $@

../bin :
	mkdir -p $@

clean :
	rm -f $(PROG) $(OBJS) $(DEPFILES) *~ gmon.out *.gcno *.gcda *.gcov \
	./utils/*.gcno ./utils/*.gcda ./utils/*.gcov \
	./arch/*/*.gcno ./arch/*/*.gcda ./arch/*/*.gcov \
	./os/*/*.gcno ./os/*/*.gcda ./os/*/*.gcov
	rm -rf ../bin/

%.o : %.cc Makefile
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -c -o $@ $<

arch/x86_64/ssse3.o : arch/x86_64/ssse3.cc Makefile
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -mssse3 -c -o $@ $<

arch/x86_64/sse41.o : arch/x86_64/sse41.cc Makefile
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -msse4.1 -c -o $@ $<

arch/x86_64/popcnt.o : arch/x86_64/popcnt.cc Makefile
	$(CXX) $(CXXFLAGS) $(DEPFLAGS) -mpopcnt -c -o $@ $<

# Pull in the auto-generated header-dependency files. The leading "-"
# silences the warning when a .d file does not yet exist.
-include $(DEPFILES)
