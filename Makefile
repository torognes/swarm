# SWARM
#
# Copyright (C) 2012 Torbjorn Rognes and Frederic Mahe
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

LIBS=-lpthread
LINKFLAGS=$(COMMON)

# Intel options
#CXX=icpc
#CXXFLAGS=-Wall -Wno-missing-declarations -fast -xSSE2 $(COMMON)

# GNU options
CXX=g++
#CXXFLAGS=-Wall -O3 -mtune=core2 -msse4.1 $(COMMON)
CXXFLAGS=-Wall -O3 -march=core2 $(COMMON)

PROG=swarm

OBJS=swarm.o db.o search8.o search16.o nw.o matrix.o util.o scan.o algo.o

DEPS=swarm.h Makefile


.SUFFIXES:.o .cc

%.o : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

%.s : %.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -c -S -o $@ $<

all : $(PROG)

swarm : $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)

clean :
	rm -f *.o *~ $(PROG) gmon.out
