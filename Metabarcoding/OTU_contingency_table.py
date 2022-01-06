#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Read all fasta files and build a sorted OTU contingency
    table. Usage: python OTU_contingency_table.py [input files]
"""

from __future__ import print_function

__author__ = "Frédéric Mahé <frederic.mahe@cirad.fr> modified by Benoit Perez-Lamarque <benoit.perez@ens.fr>"
__date__ = "2020/07/18"
__version__ = "$Revision: 5.0"

import os
import re
import sys
import operator
import numpy as np

#*****************************************************************************#
#                                                                             #
#                                  Functions                                  #
#                                                                             #
#*****************************************************************************#


def representatives_parse():
    """
    Get seed sequences.
    """
    separator = ";size="
    representatives_file = sys.argv[1]
    representatives = dict()
    with open(representatives_file, "r") as representatives_file:
        for line in representatives_file:
            if line.startswith(">"):
                amplicon = line.strip(">;\n").split(separator)[0]
            else:
                representatives[amplicon] = line.strip()

    return representatives


def stats_parse():
    """
    Map OTU seeds and stats.
    """
    separator = "\t"
    stats_file = sys.argv[2]
    stats = dict()
    with open(stats_file, "r") as stats_file:
        for line in stats_file:
            cloud, mass, seed, seed_abundance = line.strip().split(separator)[0:4]
            stats[seed] = int(mass)
    # Sort OTUs by decreasing mass
    sorted_stats = sorted(stats.items(),
                          key=operator.itemgetter(1, 0))
    sorted_stats.reverse()

    return stats, sorted_stats


def swarms_parse():
    """
    Map OTUs.
    """
    separator = "_[0-9]+|;size=[0-9]+;?| "  # parsing of abundance annotations
    swarms_file = sys.argv[3]
    swarms = dict()
    with open(swarms_file, "r") as swarms_file:
        for line in swarms_file:
            line = line.strip()
            amplicons = re.split(separator, line)[0::2]
            seed = amplicons[0]
            amplicons = [string for string in amplicons if string != '']
            #amplicons[0] = amplicons[0].strip("OTU")
            swarms[seed] = [np.unique(amplicons)]
            
    return swarms


def uchime_parse():
    """
    Map OTU's chimera status.
    """
    separator = " "
    uchime_file = sys.argv[4]
    uchime = dict()
    with open(uchime_file, "r") as uchime_file:
        for line in uchime_file:
            OTU = line.strip().split("\t")
            try:
                seed = OTU[1].split(";")[0]
            except IndexError:  # deal with partial line (missing seed)
                continue
            try:
                status = OTU[17]
            except IndexError:  # deal with unfinished chimera detection runs
                status = "NA"
            uchime[seed] = status

    return uchime



def stampa_parse():
    """
    Map amplicon ids and taxonomic assignments.
    """
    separator = "\t"
    stampa_file = sys.argv[5]
    stampa = dict()
    with open(stampa_file, "r") as stampa_file:
        for line in stampa_file:
            amplicon, identity, taxonomy = line.strip().split(separator)
            amplicon, abundance = amplicon.split(";size=")
            stampa[amplicon] = (identity, taxonomy)

    return stampa


def fasta_parse():
    """
    Map amplicon ids, abundances and samples.
    """
    separator = ";size="
    fasta_files = sys.argv[6:]
    samples = dict()
    amplicons2samples = dict()
    for fasta_file in fasta_files:
        sample = os.path.basename(fasta_file)
        sample = os.path.splitext(sample)[0]
        samples[sample] = samples.get(sample, 0) + 1
        with open(fasta_file, "r") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    amplicon, abundance = line.strip(">;\n").split(separator)
                    abundance = int(abundance)
                    if amplicon not in amplicons2samples:
                        amplicons2samples[amplicon] = {sample: abundance}
                    else:
                        amplicons2samples[amplicon][sample] = amplicons2samples[amplicon].get(sample, 0) + abundance
    # deal with duplicated samples
    duplicates = [sample for sample in samples if samples[sample] > 1]
    if duplicates:
        print("Warning: some samples are duplicated", file=sys.stderr)
        print("\n".join(duplicates), file=sys.stderr)
    samples = sorted(samples.keys())
    return amplicons2samples, samples


def print_table(representatives, stats, sorted_stats,
                swarms, uchime, amplicons2samples,
                samples, stampa):
    """
    Export results.
    """
    # Print table header
    print("OTU", "abundance",
          "amplicon", "length",
          "chimera", "spread",
          #"sequence",  # sequences makes the OTU table too heavy
          "identity", "taxonomy",
          "\t".join(samples),
          sep="\t", file=sys.stdout)

	# Print table content
    i = 1
    for seed, abundance in sorted_stats:
        sequence = representatives[seed]
        occurrences = dict([(sample, 0) for sample in samples])
        for amplicons in swarms[seed]:
            for amplicon in amplicons:
                for sample in samples:
                    occurrences[sample] += amplicons2samples[amplicon].get(sample, 0)
        spread = len([occurrences[sample] for sample in samples if occurrences[sample] > 0])


        if seed in uchime:
            chimera_status = uchime[seed]
        else:
            chimera_status = "NA"


        if seed in stampa:
            identity, taxonomy = stampa[seed]
        else:
            identity, taxonomy = "NA", "NA"

        print(i, abundance,
              seed, len(sequence),
              chimera_status, spread,
              identity, taxonomy,
              "\t".join([str(occurrences[sample]) for sample in samples]),
              sep="\t", file=sys.stdout)
        i += 1

    return


def main():
    """
    Read all fasta files and build a sorted OTU contingency table.
    """
    representatives = representatives_parse()

    stats, sorted_stats = stats_parse()

    swarms = swarms_parse()

    uchime = uchime_parse()

    stampa = stampa_parse()
    
    amplicons2samples, samples = fasta_parse()

    print_table(representatives, stats, sorted_stats, swarms,
                uchime, amplicons2samples, samples,
                stampa)

    return


#*****************************************************************************#
#                                                                             #
#                                     Body                                    #
#                                                                             #
#*****************************************************************************#

if __name__ == '__main__':

    main()

sys.exit(0)
