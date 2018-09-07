# SWARM
#
# Copyright (C) 2012-2017 Torbjorn Rognes and Frederic Mahe
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
#COMMON=-pg -g -flto -O3 -march=native -msse2
COMMON=-g -flto -O3 -march=native -msse2 -std=c++11

LIBS=-lpthread
LINKFLAGS=$(COMMON)

CXX=g++
WARNINGS=-Wall -Wsign-compare -Wextra -pedantic -Wno-long-long
CXXFLAGS=$(COMMON) $(WARNINGS) -Icityhash

PROG=swarm

OBJS=swarm.o db.o search8.o search16.o nw.o matrix.o util.o scan.o \
	algo.o algod1.o qgram.o ssse3.o derep.o arch.o cityhash/city.o \
	zobrist.o

DEPS=Makefile swarm.h bitmap.h bloom.h cityhash/config.h cityhash/city.h \
	threads.h zobrist.h

all : $(PROG)

swarm : $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -a swarm ../bin

clean :
	rm -rf swarm *.o *~ ../bin/ gmon.out cityhash/*.o \
	../man/*~ ../*~

ssse3.o : ssse3.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mssse3 -c -o $@ $<
