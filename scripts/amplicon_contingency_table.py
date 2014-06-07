#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read all fasta files and build a sorted amplicon contingency
    table. Usage: python amplicon_contingency_table.py samples_*.fas
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2014/06/03"
__version__ = "$Revision: 1.0"

import sys
import operator

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def fasta_parse():
    """
    Map amplicon ids, abundances and samples
    """
    fasta_files = sys.argv[1:]
    all_amplicons = dict()
    samples = list()
    amplicons2samples = dict()
    for fasta_file in fasta_files:
        sample = fasta_file.split(".")[0]
        samples.append(sample)
        with open(fasta_file, "rU") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    amplicon, abundance = line.strip(">\n").split("_")
                    abundance = int(abundance)
                    if amplicon not in amplicons2samples:
                        amplicons2samples[amplicon] = {sample: abundance}
                    else:
                        amplicons2samples[amplicon][sample] = abundance
                    if amplicon in all_amplicons:
                        all_amplicons[amplicon] += abundance
                    else:
                        all_amplicons[amplicon] = abundance

    return all_amplicons, amplicons2samples, samples


def main():
    """
    Read all fasta files and build a sorted amplicon contingency table
    """
    # Parse command line
    all_amplicons, amplicons2samples, samples = fasta_parse()

    # Sort amplicons by decreasing abundance (and by amplicon name)
    sorted_all_amplicons = sorted(all_amplicons.iteritems(),
                                  key=operator.itemgetter(1, 0))
    sorted_all_amplicons.reverse()

    # Print table header
    print("amplicon", "\t".join(samples), "total", sep="\t", file=sys.stdout)

    # Print table content
    for amplicon, abundance in sorted_all_amplicons:
        abundances = [amplicons2samples[amplicon].get(sample, 0)
                      for sample in samples]
        total = sum(abundances)
        abundances = [str(i) for i in abundances]
        # Sanity check
        if total == abundance:
            print(amplicon, "\t".join(abundances), total, sep="\t",
                  file=sys.stdout)
        else:
            print("Abundance sum is not correct for this amplicon",
                  amplicon, str(abundance), str(total), file=sys.stderr)
            sys.exit(-1)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
