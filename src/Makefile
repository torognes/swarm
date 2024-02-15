# SWARM
#
# Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe
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
WARNINGS = -Wall -Wextra -Wpedantic
CLANG = clang
DEBUG_WARNINGS = -Wcast-align -Wcast-qual -Wconversion -Wdate-time \
	-Wdouble-promotion -Wduplicated-branches -Wduplicated-cond -Wfloat-equal \
	-Wformat=2 -Wformat-overflow -Wlogical-op -Wnon-virtual-dtor -Wnull-dereference \
	-Wold-style-cast -Woverloaded-virtual -Wshadow -Wsign-conversion -Wuninitialized \
        -Wunsafe-loop-optimizations -Wunused -Wunused-macros -Wuseless-cast -Wvla
DEBUG_SANITIZER = -fsanitize=undefined,address -fno-sanitize=alignment -fno-omit-frame-pointer

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

# default behaviour is to compile for release (max optimizations)
else
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


# Identify Machine
ifeq ($(CXX), aarch64-linux-gnu-g++)
	MACHINE = aarch64
else ifeq ($(CXX), powerpc64le-linux-gnu-g++)
	MACHINE = ppc64le
else ifeq ($(CXX), x86_64-linux-gnu-g++)
	MACHINE = x86_64
else ifeq ($(CXX), x86_64-w64-mingw32-g++)
	MACHINE = x86_64
else
	MACHINE = $(shell uname -m)
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
ifeq ($(CXX), x86_64-w64-mingw32-g++)
	LIBS += -lpsapi
	LINKOPT += -static
	BIN = swarm.exe
else
	BIN = swarm
endif


LINKFLAGS = $(COMMON) $(LINKOPT)

CXXFLAGS = $(COMMON) $(WARNINGS)

PROG = swarm

OBJS = algo.o algod1.o arch.o bloomflex.o bloompat.o db.o derep.o \
	hashtable.o nw.o qgram.o scan.o search16.o search8.o \
	swarm.o util.o utils/cigar.o utils/fatal.o utils/hashtable_size.o \
	utils/input_output.o \
	utils/intrinsics_to_functions_aarch64.o utils/intrinsics_to_functions_ppc.o \
	utils/intrinsics_to_functions_x86_64.o \
	utils/nt_codec.o utils/progress.o \
	variants.o x86_cpu_features.o zobrist.o $(EXTRAOBJ)

DEPS = Makefile arch.h algo.h algod1.h bloomflex.h bloompat.h db.h derep.h hashtable.h \
	nw.h qgram.h scan.h search16.h search8.h swarm.h \
	utils/alignment_parameters.h \
	util.h utils/backtrack.h utils/cigar.h utils/fatal.h utils/gcd.h utils/hashtable_size.h \
	utils/input_output.h \
	utils/intrinsics_to_functions_aarch64.h utils/intrinsics_to_functions_ppc.h \
	utils/intrinsics_to_functions_x86_64.h \
	utils/nt_codec.h utils/opt_boundary.h utils/progress.h utils/pseudo_rng.h \
	utils/qgram_array.h utils/qgram_threadinfo.h utils/queryinfo.h utils/search_data.h \
	utils/score_matrix.h utils/threads.h utils/x86_cpu_feature_popcnt.h variants.h x86_cpu_features.h \
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
