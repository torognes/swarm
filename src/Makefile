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
#COMMON=-pg -g
COMMON=-g

COMPILEOPT=-Wall -Wsign-compare -O3 -msse2 -mtune=core2 -Icityhash

LIBS=-lpthread
LINKFLAGS=$(COMMON)

CXX=g++
CXXFLAGS=$(COMPILEOPT) $(COMMON)

PROG=swarm

OBJS=swarm.o db.o search8.o search16.o nw.o matrix.o util.o scan.o \
	algo.o algod1.o qgram.o ssse3.o derep.o arch.o cityhash/city.o

DEPS=Makefile swarm.h bitmap.h bloom.h cityhash/config.h cityhash/city.h \
	threads.h

.SUFFIXES:.o .cc

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

all : $(PROG)

swarm : $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -a swarm ../bin

clean :
	rm -rf swarm *.o *~ ../bin/ gmon.out cityhash/*.o ../man/*~ ../*~

ssse3.o : ssse3.cc $(DEPS)
	$(CXX) -mssse3 $(CXXFLAGS) -c -o $@ $<
