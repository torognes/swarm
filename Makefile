# SWARM
#
# Copyright (C) 2012-2026 Torbjorn Rognes and Frederic Mahe
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

PROG := bin/swarm
MAN  := man/swarm.1
SRC  := src
BASH_COMPLETION := completion/swarm.bash
ZSH_COMPLETION  := completion/_swarm

PREFIX ?= /usr/local
exec_prefix := $(PREFIX)
datarootdir := $(PREFIX)/share
bindir      := $(exec_prefix)/bin
mandir      := $(datarootdir)/man
man1dir     := $(mandir)/man1
bashcompdir ?= $(datarootdir)/bash-completion/completions
zshcompdir  ?= $(datarootdir)/zsh/site-functions

INSTALL         ?= /usr/bin/install
INSTALL_PROGRAM ?= $(INSTALL) -m 0755
INSTALL_DATA    ?= $(INSTALL) -m 0644
MKDIR_P         ?= $(INSTALL) -d
RM              ?= rm -f

.PHONY: all swarm install install-completion uninstall clean distclean $(PROG)

all: swarm

swarm: $(PROG)

$(PROG):
	$(MAKE) -C $(SRC)

install: $(PROG) $(MAN) install-completion
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(PROG) $(DESTDIR)$(bindir)
	$(MKDIR_P) $(DESTDIR)$(man1dir)
	$(INSTALL_DATA) $(MAN) $(DESTDIR)$(man1dir)

install-completion: $(BASH_COMPLETION) $(ZSH_COMPLETION)
	$(MKDIR_P) $(DESTDIR)$(bashcompdir)
	$(INSTALL_DATA) $(BASH_COMPLETION) $(DESTDIR)$(bashcompdir)/swarm
	$(MKDIR_P) $(DESTDIR)$(zshcompdir)
	$(INSTALL_DATA) $(ZSH_COMPLETION) $(DESTDIR)$(zshcompdir)/_swarm

uninstall:
	$(RM) $(DESTDIR)$(bindir)/$(notdir $(PROG))
	$(RM) $(DESTDIR)$(man1dir)/$(notdir $(MAN))
	$(RM) $(DESTDIR)$(bashcompdir)/swarm
	$(RM) $(DESTDIR)$(zshcompdir)/_swarm

clean:
	$(MAKE) -C $(SRC) clean

distclean: clean
	$(RM) *~ $(SRC)/*~ $(SRC)/utils/*~ man/*~
