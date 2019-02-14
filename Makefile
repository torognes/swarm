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

PROG=bin/swarm
MAN=man/swarm.1

swarm : $(PROG)

$(PROG) :
	make -C src swarm

install : $(PROG) $(MAN)
	/usr/bin/install -c $(PROG) '/usr/local/bin'
	/usr/bin/install -c $(MAN) '/usr/local/share/man/man1'

clean :
	make -C src clean
