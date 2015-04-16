#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Visualize the internal structure of a swarm (color vertices by
    abundance). Amplicons grafted with the fastidious option will be
    discarded and will not be visualized. Requires the module igraph
    and python 2.7+.
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <mahe@rhrk.uni-kl.fr>"
__date__ = "2015/04/16"
__version__ = "$Revision: 2.0"

import sys
import os.path
from igraph import Graph, plot
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
    desc = """Visualize the internal structure of a given OTU."""

    parser = OptionParser(usage="usage: %prog -s FILE -i FILE -o INT",
                          description=desc,
                          version="%prog version 2.0")

    parser.add_option("-s", "--swarms",
                      metavar="<FILENAME>",
                      action="store",
                      dest="swarms",
                      help="<FILENAME> contains swarm's results")

    parser.add_option("-i", "--internal_structure",
                      metavar="<FILENAME>",
                      action="store",
                      dest="internal_structure",
                      help="<FILENAME> contains OTUs' internal structure")

    parser.add_option("-o", "--OTU",
                      metavar="<INTEGER>",
                      action="store",
                      type="int",
                      dest="OTU",
                      default=1,
                      help="Select the nth OTU (first by default)")

    (options, args) = parser.parse_args()

    return options.swarms, options.internal_structure, options.OTU


def parse_files(swarms, internal_structure, OTU):
    """
    """
    # List amplicon ids and abundances
    amplicons = list()
    with open(swarms, "rU") as swarms:
        for i, swarm in enumerate(swarms):
            if i == OTU - 1:
                amplicons = [tuple(item.split("_"))
                             for item in swarm.strip().split(" ")]
                break

    # Convert amplicon names to amplicon indexes (igraph uses
    # numerical ids, not names). Requires python 2.7+)
    amplicon_index = {amplicon[0]: i for (i, amplicon) in enumerate(amplicons)}
    amplicon_connected = {amplicon[0]: False for amplicon in amplicons}

    # List pairwise relations
    relations = list()
    with open(internal_structure, "rU") as internal_structure:
        for line in internal_structure:
            # Get the first four elements of the line
            ampliconA, ampliconB, d, OTU_number = line.strip().split("\t")[0:4]
            if int(OTU_number) == OTU:
                amplicon_connected[ampliconA] = True
                amplicon_connected[ampliconB] = True
                relations.append((amplicon_index[ampliconA],
                                  amplicon_index[ampliconB]))

    # Discard amplicons grafted with the fastidious option
    amplicons = [amplicon for amplicon in amplicons
                 if amplicon_connected[amplicon[0]]]

    # Deal with errors
    if not amplicons or not relations:
        print("Error: OTU does not exists or contains only one element.",
              file=sys.stderr)
        sys.exit(-1)

    return amplicons, relations


def build_graph(amplicons, relations):
    """
    Convert pairwise relations into a graph structure.
    """
    # Create vertices (= number of unique amplicons)
    g = Graph(len(amplicons))
    g.add_edges(relations)

    amplicon_ids = [amplicon[0] for amplicon in amplicons]
    abundances = [int(amplicon[1]) for amplicon in amplicons]
    minimum, maximum = min(abundances), max(abundances)

    # Determine canvas size
    if len(abundances) < 500:
        bbox = (1920, 1080)
    elif len(abundances) > 4000:
        bbox = (5760, 3240)
    else:
        bbox = (3840, 2160)

    # Compute node attributes
    node_colors = list()
    node_sizes = list()
    node_labels = list()
    for abundance in abundances:
        # Color is coded by a 3-tuple of float values (0.0 to 1.0)
        # Start from a max color in rgb(red, green, blue)
        max_color = (176, 196, 222)  # light steel blue
        color = [1.0 * (c + (255 - c) / abundance) / 255 for c in max_color]
        node_colors.append(color)
        node_size = 30 + (abundance * 70 / maximum)
        node_sizes.append(node_size)
        # Label nodes with an abundance greater than 10
        if abundance >= 10 or abundance == maximum:
            node_labels.append(str(abundance))
        else:
            node_labels.append("")  # Doesn't work with "None"

    g.vs["name"] = amplicon_ids
    g.vs["abundance"] = abundances
    g.vs["label"] = node_labels
    g.vs["color"] = node_colors
    g.vs["vertex_size"] = node_sizes

    return g, bbox


def main():
    """
    Backbone of the script
    """
    # Parse command line options.
    swarms, internal_structure, OTU = option_parse()

    # Output file
    basename = os.path.splitext(swarms)[0]
    output_pdf = basename + "_OTU_" + str(OTU) + ".pdf"
    # output_svg = basename + "_OTU_" + str(OTU) + ".svg"
    output_graphml = basename + "_OTU_" + str(OTU) + ".graphml"

    # Collect data
    amplicons, relations = parse_files(swarms, internal_structure, OTU)

    # Build the graph with igraph
    g, bbox = build_graph(amplicons, relations)

    # Select the attributes to display
    visual_style = dict()
    visual_style["vertex_size"] = g.vs["vertex_size"]
    visual_style["vertex_label"] = g.vs["label"]
    visual_style["vertex_color"] = g.vs["color"]
    visual_style["vertex_label_dist"] = 0
    visual_style["layout"] = g.layout("fr")  # alternatives: "kk", "drl"
    visual_style["bbox"] = bbox
    visual_style["margin"] = 20

    # Output the graph
    g.write_graphml(output_graphml)
    plot(g, output_pdf, **visual_style)
    # plot(g, output_svg, **visual_style)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
