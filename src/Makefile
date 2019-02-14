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

LIBS=-lpthread
LINKFLAGS=$(COMMON)

CXX=g++
WARNINGS=-Wall -Wsign-compare -Wextra -pedantic -Wno-long-long
CXXFLAGS=$(COMMON) $(WARNINGS) -O3 -msse2 -mtune=core2 -Icityhash

PROG=swarm
MAN=../man/swarm.1

SRCS=swarm.cc db.cc search8.cc search16.cc nw.cc matrix.cc util.cc scan.cc \
	algo.cc algod1.cc qgram.cc ssse3.cc derep.cc arch.cc cityhash/city.cc

OBJS=$(SRCS:.cc=.o)

DEPS=Makefile swarm.h bitmap.h bloom.h cityhash/config.h cityhash/city.h \
	threads.h

all : $(PROG)

swarm : $(OBJS)
	$(CXX) $(LINKFLAGS) -o $@ $(OBJS) $(LIBS)
	mkdir -p ../bin
	cp -a swarm ../bin
	gzip -kf $(MAN)

ssse3.o : ssse3.cc $(DEPS)
	$(CXX) $(CXXFLAGS) -mssse3 -c -o $@ $<

install : $(PROG) $(MAN)
	/usr/bin/install -c $(PROG) '/usr/local/bin'
	/usr/bin/install -c $(MAN:.1=.1.gz) '/usr/local/share/man/man1'

clean :
	rm -rf swarm *.o *~ ../bin/ gmon.out cityhash/*.o ../man/*~ ../man/*.gz ../*~
