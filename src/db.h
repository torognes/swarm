/*
    SWARM

    Copyright (C) 2012-2021 Torbjorn Rognes and Frederic Mahe

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
    Department of Informatics, University of Oslo,
    PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include <cstdio>
#include <cstdint>

void db_read(const char * filename, struct Parameters const & p);

auto db_getsequencecount() -> unsigned int;

auto db_getnucleotidecount() -> uint64_t;

auto db_getlongestheader() -> unsigned int;

auto db_getlongestsequence() -> unsigned int;

auto db_getseqinfo(uint64_t seqno) -> struct seqinfo_s *;

auto db_getsequence(uint64_t seqno) -> char *;

auto db_getsequencelen(uint64_t seqno) -> unsigned int;

auto db_gethash(uint64_t seqno) -> uint64_t;

void db_getsequenceandlength(uint64_t seqno,
                             char ** address,
                             unsigned int * length);

auto db_getheader(uint64_t seqno) -> char *;

auto db_getheaderlen(uint64_t seqno) -> unsigned int;

auto db_getabundance(uint64_t seqno) -> uint64_t;

void db_showall();

void db_free();

void db_putseq(int64_t seqno);

void db_qgrams_init();

void db_qgrams_done();

void db_fprintseq(std::FILE * fp, unsigned int a, unsigned int width);

void fprint_id(std::FILE * stream,
               uint64_t x,
               bool opt_usearch_abundance,
               int64_t opt_append_abundance);

void fprint_id_noabundance(std::FILE * stream,
                           uint64_t x,
                           bool opt_usearch_abundance);

void fprint_id_with_new_abundance(std::FILE * stream,
                                  uint64_t seqno,
                                  uint64_t abundance,
                                  bool opt_usearch_abundance);
