# SWARM
#
# Copyright (C) 2012-2023 Torbjorn Rognes and Frederic Mahe
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

COMMON = -g
EXTRAOBJ =
LINKOPT =
LIBS = -lpthread
WARNINGS = -Wall -Wextra
CLANG = clang

# Run "make RELEASE=1" to compile for release
ifdef RELEASE
	COMMON += -O3 -DNDEBUG
	COMPILER_VERSION := $(shell $(CXX) --version)
	ifeq ($(CLANG), $(findstring $(CLANG), $(COMPILER_VERSION)))
		COMMON += -flto=full
	else
		GCCVERSIONGTEQ12 := $(shell expr `gcc -dumpversion | cut -f 1 -d "."` \>= 12)
		ifeq "$(GCCVERSIONGTEQ12)" "1"
			COMMON += -flto=auto
		else
			COMMON += -flto
		endif
	endif
endif

# Run "make PROFILE=1" to compile for profiling
ifdef PROFILE
	COMMON += -pg
endif

# Run "make DEBUG=1" to compile for debugging
ifdef DEBUG
	COMMON += -Og -DDEBUG
endif

# Run "make COVERAGE=1" to compile for coverage tests
ifdef COVERAGE
	COMMON += -fprofile-arcs -ftest-coverage -O0
	LIBS += -lgcov
	LINKFLAGS += --coverage
endif

# Machine specific
ifeq ($(CXX), aarch64-linux-gnu-g++)
	MACHINE = aarch64
else ifeq ($(CXX), powerpc64le-linux-gnu-g++)
	MACHINE = ppc64le
else
	MACHINE = $(shell uname -m)
endif

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

# OS specific
ifeq ($(CXX), x86_64-w64-mingw32-g++)
	LIBS += -lpsapi
	LINKOPT += -static
	BIN = swarm.exe
else
	WARNINGS += -pedantic
	BIN = swarm
endif

LINKFLAGS = $(COMMON) $(LINKOPT)

CXXFLAGS = $(COMMON) $(WARNINGS)

PROG = swarm

OBJS = algo.o algod1.o arch.o bloomflex.o bloompat.o db.o derep.o \
	hashtable.o nw.o qgram.o scan.o search16.o search8.o \
	swarm.o util.o utils/fatal.o utils/hashtable_size.o variants.o \
	x86_cpu_features.o zobrist.o $(EXTRAOBJ)

DEPS = Makefile arch.h algo.h algod1.h bloomflex.h bloompat.h db.h derep.h hashtable.h \
	nw.h qgram.h scan.h search16.h search8.h swarm.h \
	util.h utils/fatal.h utils/gcd.h utils/hashtable_size.h utils/nt_codec.h utils/pseudo_rng.h \
	utils/qgram_array.h utils/score_matrix.h utils/threads.h variants.h x86_cpu_features.h \
	zobrist.h


all : $(PROG)

swarm : $(OBJS) $(DEPS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -a $(BIN) ../bin

clean :
	rm -f swarm *.o ./utils/*.o *~ gmon.out *.gcno *.gcda *.gcov \
	./utils/*.gcno ./utils/*.gcda ./utils/*.gcov ../bin/swarm

.o : .cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ssse3.o : ssse3.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mssse3 -c -o $@ $<

sse41.o : sse41.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -msse4.1 -c -o $@ $<

popcnt.o : popcnt.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mpopcnt -c -o $@ $<
