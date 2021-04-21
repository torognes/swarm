# SWARM
#
# Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe
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

# Run "make RELEASE=1" to compile for release
ifdef RELEASE
	COMMON += -DNDEBUG
endif

# Run "make PROFILE=1" to compile for profiling
ifdef PROFILE
	COMMON += -pg
endif

# Run "make COVERAGE=1" to compile for coverage tests
ifdef COVERAGE
	COMMON += -fprofile-arcs -ftest-coverage -O0
	LIBS += -lgcov
	LINKFLAGS += --coverage
else
	COMMON += -flto -O3
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
	COMMON += -mcpu=power8 -std=gnu++11
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

OBJS = swarm.o db.o search8.o search16.o nw.o matrix.o util.o scan.o \
	algo.o algod1.o qgram.o derep.o arch.o \
	zobrist.o bloompat.o bloomflex.o variants.o hashtable.o \
	x86_cpu_features.o \
	$(EXTRAOBJ)

DEPS = Makefile swarm.h \
	threads.h zobrist.h bloompat.h bloomflex.h variants.h hashtable.h \
	x86_cpu_features.h matrix.h db.h


all : $(PROG)

swarm : $(OBJS) $(DEPS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -a $(BIN) ../bin

clean :
	rm -f swarm *.o *~ gmon.out *.gcno *.gcda *.gcov ../bin/swarm

.o : .cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ssse3.o : ssse3.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mssse3 -c -o $@ $<

sse41.o : sse41.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -msse4.1 -c -o $@ $<

popcnt.o : popcnt.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mpopcnt -c -o $@ $<
