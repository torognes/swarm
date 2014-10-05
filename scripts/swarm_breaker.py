#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Detect and break chains of amplicons in a swarm.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2014/01/30"
__version__ = "$Revision: 1.1"

import os
import sys
import tempfile
import itertools
import subprocess
from operator import itemgetter
from optparse import OptionParser

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def option_parse():
    """
    Parse arguments from command line.
    """
    desc = """Detect and break chains of amplicons in a swarm. That
    script will search for the swarm binary in the directories listed
    under $PATH. If swarm is installed at a different location, please
    modify the corresponding line in the function run_swarm, or use
    the -b option to specify the path to swarm."""

    parser = OptionParser(usage="usage: %prog -f filename -s filename",
                          description=desc,
                          version="%prog version 1.1")

    parser.add_option("-b", "--binary",
                      metavar="<BINARY>",
                      action="store",
                      default="swarm",
                      dest="binary",
                      help="swarm binary location. Default is /usr/bin/swarm")

    parser.add_option("-f", "--fasta_file",
                      metavar="<FILENAME>",
                      action="store",
                      dest="fasta_file",
                      help="set <FILENAME> as fasta file.")

    parser.add_option("-s", "--swarm_file",
                      metavar="<FILENAME>",
                      action="store",
                      dest="swarm_file",
                      help="set <FILENAME> as swarm file.")

    parser.add_option("-d", "--differences",
                      metavar="<THRESHOLD>",
                      action="store",
                      type="int",
                      default=1,
                      dest="threshold",
                      help="set local clustering <THRESHOLD>. Default is 1")

    parser.add_option("-t", "--threads",
                      metavar="<THREADS>",
                      action="store",
                      type="int",
                      default=1,
                      dest="threads",
                      help="set the number of <THREADS>. Default is 1")

    parser.add_option("-z", "--usearch_abundance",
                      action="store_true",
                      dest="usearch_style",
                      default=False,
                      help="abundance annotation in usearch style")

    (options, args) = parser.parse_args()
    return options.binary, options.fasta_file, options.swarm_file, \
        options.threshold, options.threads, options.usearch_style


def check_files(paths):
    """
    Check input file and swarm binary status
    """
    status = dict()
    for path in paths:
        status[path] = {"exist": True, "read": True, "execute": True}
        if not os.access(path, os.F_OK):
            status[path]["exist"] = False
        if not os.access(path, os.R_OK):
            status[path]["read"] = False
        if not os.access(path, os.X_OK):
            status[path]["execute"] = False
    # Is the path indicated by the user correct?
    if paths[0] != "swarm":
        if status[paths[0]]["exist"] is False:
            print("Error: ",
                  "Cannot find the swarm binary at ",
                  paths[0], "\n",
                  "Use -b /path/to/swarm to indicate the correct path.",
                  sep="", file=sys.stderr)
            sys.exit(-1)
        if status[paths[0]]["execute"] is False:
            print("Error: ",
                  "Cannot execute the swarm binary at ",
                  paths[0], "\n",
                  "please check file permissions.",
                  sep="", file=sys.stderr)
            sys.exit(-1)
    # Are the other files readable?
    for path in paths[1:]:
        if status[path]["exist"] is False:
            print("ERROR", "Cannot find the file", path,
                  "\nExit", file=sys.stderr)
            sys.exit(-1)
        if status[path]["read"] is False:
            print("ERROR", "Cannot read the file", path,
                  "\nExit", file=sys.stderr)
            sys.exit(-1)
    return status


def fasta_parse(fasta_file, usearch_style):
    """
    Parse the fasta file (linearized or not). List amplicon ids,
    abundances and sequences, make a dictionary
    """
    with open(fasta_file, "rU") as fasta_file:
        # What is the amplicon-abundance separator?
        separator = "_"
        if usearch_style:
            separator = ";size="

        all_amplicons = dict()
        sequence = list()
        amplicon = ""
        abundance = 0
        # Parse the fasta file
        for line in fasta_file:
            line = line.strip()
            if not line:  # skip empty lines
                continue
            if line.startswith(">"):
                if sequence:  # then finalize the previous fasta entry
                    full_sequence = "".join(sequence)
                    all_amplicons[amplicon] = (int(abundance), full_sequence)
                    sequence = list()
                # store the header and the abundance value of the new entry
                amplicon, abundance = line.lstrip(">").rsplit(separator, 1)
            else:
                sequence.append(line)
        else:  # deal with the last entry
            full_sequence = "".join(sequence)
            all_amplicons[amplicon] = (int(abundance), full_sequence)
        return all_amplicons


def swarm_parse(swarm_file):
    """
    List amplicons contained in each swarms, sort by decreasing
    abundance. Sort the list of swarms by decreasing mass and
    decreasing size.
    """
    with open(swarm_file, "rU") as swarm_file:
        swarms = list()
        for line in swarm_file:
            amplicons = [(amplicon.rsplit("_", 1)[0], int(amplicon.rsplit("_", 1)[1]))
                         for amplicon in line.strip().split(" ")]
            # Sort amplicons by decreasing abundance and alphabetical order
            amplicons.sort(key=itemgetter(1, 0), reverse=True)
            top_amplicon, top_abundance = amplicons[0]
            swarm_size = len(amplicons)
            swarm_mass = sum([amplicon[1] for amplicon in amplicons])
            swarms.append([top_amplicon, swarm_mass, swarm_size,
                           top_abundance, amplicons])
        # Sort swarms on mass, size and seed name
        swarms.sort(key=itemgetter(1, 2, 0), reverse=True)
        return swarms


def run_swarm(binary, all_amplicons, swarm, threshold, threads):
    """
    Write temporary fasta files, run swarm and collect the graph data
    """
    swarm_command = [binary, "-b", "-d", str(threshold), "-t", str(threads)]
    if threshold == 1:
        swarm_command.insert(1, "-a")  # Use the fast algorithm if d = 1
    with open(os.devnull, "w") as devnull:
        with tempfile.SpooledTemporaryFile() as tmp_fasta_file:
            with tempfile.SpooledTemporaryFile() as tmp_swarm_results:
                for amplicon, abundance in swarm:
                    sequence = all_amplicons[amplicon][1]
                    print(">", amplicon, "_", str(abundance), "\n", sequence,
                          sep="", file=tmp_fasta_file)
                tmp_fasta_file.seek(0)  # rewind to the begining of the file
                try:
                    proc = subprocess.Popen(swarm_command,
                                            stderr=tmp_swarm_results,
                                            stdout=devnull,
                                            stdin=tmp_fasta_file,
                                            close_fds=True)
                except (OSError, 2):
                    print("Error ",
                          "Cannot find the swarm binary in $PATH\n",
                          "Use -b /path/to/swarm to indicate the correct path.",
                          sep="", file=sys.stderr)
                    sys.exit(-1)
                proc.wait()  # Important to wait
                tmp_swarm_results.seek(0)  # rewind to the begining of the file
                graph_data = [line.strip().split("\t")[1:4]
                              for line in tmp_swarm_results
                              if "@@" in line]
                return graph_data


def build_graph(graph_data):
    """
    List pairwise relations in a swarm. Note that not all pairwise
    relations are stored. That's why the graph exploration must always
    start from the most abundant amplicon, and must be reiterated for
    sub-swarms after a breaking.
    """
    graph = dict()
    for line in graph_data:
        ampliconA, ampliconB, differences = line
        if ampliconA in graph:
            graph[ampliconA] += [ampliconB]
        else:
            graph[ampliconA] = [ampliconB]
    return graph


def find_path(graph, start, end, path=[]):
    """
    Recursively explore the graph and find all paths connecting two
    amplicons (http://www.python.org/doc/essays/graphs.html). As the
    graph is not complete, some pairs of amplicon cannot be linked.
    """
    path = path + [start]
    if start == end:
        return path
    if start not in graph:
        return None
    for node in graph[start]:
        if node not in path:
            newpath = find_path(graph, node, end, path)
            if newpath:
                return newpath
    return None


def graph_breaker(amplicons, graph, all_amplicons, ABUNDANT):
    """
    Find deep valleys and cut the graph
    """
    # High peaks to test (starting and ending points)
    top_amplicons = [amplicon[0] for amplicon in amplicons
                     if amplicon[1] >= ABUNDANT]
    # Ending peak is RATIO times higher than the valley
    RATIO = 50
    # Debugging
    print("## OTU ", top_amplicons[0], "\n", "# List potential bridges",
          sep="", file=sys.stderr)
    # Initialize the list of new seeds
    new_swarm_seeds = [top_amplicons[0]]
    # Break if there is no second peak
    if len(top_amplicons) < 2:
        return new_swarm_seeds, graph
    # Loop over the list of top amplicons
    pairs_of_peaks = itertools.combinations(top_amplicons, 2)
    for pair_of_peaks in pairs_of_peaks:
        start_amplicon, end_amplicon = pair_of_peaks
        path = find_path(graph, start_amplicon, end_amplicon)
        # Path can be empty if the relation have been deleted
        if path and len(path) > 1:
            abundances = [int(all_amplicons[node][0]) for node in path]
            # Find the weakest spot
            lowest = min(abundances)
            if lowest != abundances[-1]:
                # LOW VALLEY MODEL (CHANGE HERE)
                peak1 = abundances[0]
                peak2 = abundances[-1]
                if (peak2 / lowest > RATIO / 2 and peak1 / peak2 < 10) or peak2 / lowest >= RATIO:
                    # Debugging
                    print(abundances, "\tBREAK!", file=sys.stderr)
                    # Find the rightmost occurence of the lowest point
                    index = len(abundances) - (abundances[::-1].index(lowest) + 1)
                    left_amplicon = path[index-1]
                    right_amplicon = path[index]
                    # Delete the relation from the graph
                    graph[left_amplicon].remove(right_amplicon)
                    # Remove the graph entry if the relation is now empty
                    if not graph[left_amplicon]:
                        del graph[left_amplicon]
                    # Lowest point will be a new swarm seed
                    new_swarm_seeds.append(right_amplicon)
                else:
                    print(abundances, file=sys.stderr)
    return new_swarm_seeds, graph


def swarmer(graph, seed, path=[]):
    """
    Recursively explore the graph and find all amplicons linked to the
    seed
    """
    path = path + [seed]
    if seed in graph:
        for node in graph[seed]:
            path = swarmer(graph, node, path)
    return path


def swarm_breaker(binary, all_amplicons, swarms, threshold, threads):
    """
    Recursively inspect and break the newly produced swarms
    """
    # ARBITRARY PARAMETERS
    ABUNDANT = 100
    # Deal with each swarm
    for swarm in swarms:
        top_amplicon, swarm_mass, swarm_size, top_abundance, amplicons = swarm
        if swarm_size > 2 and top_abundance > ABUNDANT:
            # Run swarm to get the pairwise relationships
            graph_raw_data = run_swarm(binary, all_amplicons,
                                       amplicons, threshold, threads)
            # Build the graph of pairwise relationships
            graph = build_graph(graph_raw_data)
            new_swarm_seeds, graph = graph_breaker(amplicons, graph,
                                                   all_amplicons, ABUNDANT)
            # Explore the graph and find all amplicons linked to the seeds
            observed = 0
            new_swarms = list()
            for seed in new_swarm_seeds:
                new_swarm = swarmer(graph, seed)
                observed += len(new_swarm)
                # Give to the new swarms the same structure and
                # re-order them by decreasing abundance
                amplicons = [(amplicon, all_amplicons[amplicon][0])
                             for amplicon in new_swarm]
                amplicons.sort(key=itemgetter(1), reverse=True)
                top_amplicon, top_abundance = amplicons[0]
                swarm_size = len(amplicons)
                swarm_mass = sum([amplicon[1] for amplicon in amplicons])
                new_swarms.append([top_amplicon, swarm_mass,
                                   swarm_size, top_abundance, amplicons])
            # Deal with the new swarms (no need to treat again the
            # first swarm). There will always be at least one swarm in
            # new_swarms.
            print(" ".join(["_".join([amplicon[0], str(amplicon[1])])
                            for amplicon in new_swarms[0][4]]),
                  file=sys.stdout)
            new_swarms.pop(0)
            if new_swarms:
                # Sort the rest of the new swarms by decreasing mass
                # and size. Inject them into swarm_breaker.
                new_swarms.sort(key=itemgetter(1, 2), reverse=True)
                swarm_breaker(binary, all_amplicons,
                              new_swarms, threshold, threads)
        else:
            # Output the swarm
            print(" ".join(["_".join([amplicon[0], str(amplicon[1])])
                            for amplicon in amplicons]), file=sys.stdout)
    return None


def main():
    """
    Hypothesis: chain of amplicons happen among the most abundant
    amplicons of the swarm. The number of chains in a swarm is
    small. The abundances of each sub-swarm centroids are
    comparable. The "valleys" are deep compared to the "peaks". Swarm
    graphs are acyclical, so there is only one path joining two
    amplicons.

    Synopsis: Break bridges as you discover them. Find the weakest
    point in the chain, break on the left of that point and mark it as
    the seed of a new swarm. Repeat the process with the nth most
    abundant amplicon, until all amplicons in the arbitrary range have
    been treated.
    """
    # Parse command line options.
    binary, fasta_file, swarm_file, threshold, \
        threads, usearch_style = option_parse()
    # Check files
    check_files([binary, fasta_file, swarm_file])
    # Load all amplicon ids, abundances and sequences
    all_amplicons = fasta_parse(fasta_file, usearch_style)
    # Load the swarming data
    swarms = swarm_parse(swarm_file, usearch_style)
    # Deal with each swarm
    swarm_breaker(binary, all_amplicons, swarms, threshold, threads)


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
