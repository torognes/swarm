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
EXTRAOBJ =
LINKOPT =
LIBS =
WARNINGS = -Wall -Wextra -Wpedantic
DEBUG_WARNINGS = -Wcast-align -Wcast-qual -Wconversion -Wdate-time \
	-Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wfloat-equal \
	-Wformat=2 -Wformat-overflow -Wlogical-op -Wnon-virtual-dtor -Wnull-dereference \
	-Wold-style-cast -Woverloaded-virtual -Wshadow -Wsign-conversion -Wuninitialized \
        -Wunsafe-loop-optimizations -Wunused -Wunused-macros -Wuseless-cast -Wvla
DEBUG_SANITIZER = -fsanitize=undefined,address -fno-omit-frame-pointer


# Compiler identity and major version. Use the preprocessor to detect
# clang via the __clang__ macro: that works regardless of the binary
# name (g++, clang++, ccache wrappers, cross-compilers...). Use $(CXX)
# -dumpversion to query the actual compiler being used, not the host
# /usr/bin/gcc.
IS_CLANG := $(shell $(CXX) -x c++ -E -dM - < /dev/null 2>/dev/null | grep -c '__clang__')
CXX_VERSION_MAJOR := $(shell $(CXX) -dumpversion 2>/dev/null | cut -d. -f1)
CXX_VERSION_MAJOR := $(or $(CXX_VERSION_MAJOR),0)


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
	LIBS += -lgcov
	LINKFLAGS += --coverage
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
ifeq ($(MACHINE), x86_64)
	COMMON += -march=x86-64 -mtune=generic -std=c++11
	EXTRAOBJ += ssse3.o sse41.o popcnt.o
else ifeq ($(MACHINE), aarch64)
	COMMON += -march=armv8-a+simd -mtune=generic \
	          -flax-vector-conversions -std=c++11
else ifeq ($(MACHINE), arm64)
	COMMON += -std=c++11
else ifeq ($(MACHINE), ppc64le)
	COMMON += -mcpu=power8 -maltivec -std=gnu++11
endif


# OS specific (Windows vs. others)
ifneq (,$(findstring mingw, $(TARGET_TRIPLE)))
	LIBS += -lpsapi
	LINKOPT += -static
	BIN = swarm.exe
else
	BIN = swarm
endif


LINKFLAGS = $(COMMON) $(LINKOPT)

CXXFLAGS = $(COMMON) $(WARNINGS)

PROG = swarm

OBJS = algo.o algod1.o arch.o cli.o db.o derep.o \
	nw.o qgram.o scan.o search16.o search8.o \
	swarm.o variants.o \
	$(patsubst %.cc, %.o, $(wildcard utils/*.cc)) $(EXTRAOBJ)

DEPS = Makefile $(wildcard *.h) $(wildcard utils/*.h)

.PHONY: all swarm clean

all : $(PROG)

swarm : $(OBJS) $(DEPS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -a $(BIN) ../bin

clean :
	rm -f $(BIN) $(OBJS) *~ gmon.out *.gcno *.gcda *.gcov \
	./utils/*.gcno ./utils/*.gcda ./utils/*.gcov
	rm -rf ../bin/

.o : .cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ssse3.o : ssse3.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mssse3 -c -o $@ $<

sse41.o : sse41.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -msse4.1 -c -o $@ $<

popcnt.o : popcnt.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mpopcnt -c -o $@ $<
