# SWARM
#
# Copyright (C) 2012-2019 Torbjorn Rognes and Frederic Mahe
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

# Profiling options
# PROFILING=-pg
PROFILING=

# Machine specific
MACHINE=$(shell uname -m)
ifeq ($(MACHINE), x86_64)
	ARCHOPT = -march=x86-64 -mtune=generic -std=c++11
	EXTRAOBJ = ssse3.o
else ifeq ($(MACHINE), aarch64)
	ARCHOPT = -march=armv8-a+simd -mtune=generic \
	          -flax-vector-conversions -std=c++11
	EXTRAOBJ =
else ifeq ($(MACHINE), ppc64le)
	ARCHOPT = -mcpu=power8 -std=gnu++11
	EXTRAOBJ =
endif

# OS specific
ifeq ($(CXX), x86_64-w64-mingw32-g++)
	LIBS = -lpthread -lpsapi
	WARNOPT =
	LINKOPT = -static
else
	LIBS = -lpthread
	WARNOPT = -pedantic
	LINKOPT =
endif

WARNINGS = -Wall -Wsign-compare -Wextra -Wno-long-long $(WARNOPT) -Wshorten-64-to-32 -Weverything -Wno-c++98-compat -Wno-c++98-compat-pedantic -Wno-padded

COMMON=$(PROFILING) -g -flto -O3 $(ARCHOPT)

LINKFLAGS=$(COMMON) $(LINKOPT)

CXXFLAGS=$(COMMON) $(WARNINGS)

PROG=swarm

OBJS=swarm.o db.o search8.o search16.o nw.o matrix.o util.o scan.o \
	algo.o algod1.o qgram.o derep.o arch.o city.o \
	zobrist.o bloompat.o bloomflex.o variants.o hashtable.o \
	$(EXTRAOBJ)

DEPS=Makefile swarm.h city.h citycrc.h \
	threads.h zobrist.h bloompat.h bloomflex.h variants.h hashtable.h

all : $(PROG)

swarm : $(OBJS) $(DEPS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -a swarm ../bin

clean :
	rm -rf swarm *.o *~ gmon.out

.o : .cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

ssse3.o : ssse3.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mssse3 -c -o $@ $<
