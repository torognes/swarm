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

# Makefile for Swarm PDF manual from Swarm man file

all : swarm_manual.pdf

swarm_manual.pdf : swarm.1
	sed -e 's/\\-/-/g' $< | \
	iconv -f UTF-8 -t ISO-8859-1 | \
	groff -t -m mandoc -T ps -P -pa4 | ps2pdf - $@

clean : swarm_manual.pdf
	rm swarm_manual.pdf
