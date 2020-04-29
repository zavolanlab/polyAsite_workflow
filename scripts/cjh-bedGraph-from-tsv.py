__date__ = "2019-07-24"
__author__ = "Christina J. Herrmann"
__email__ = "christina.herrmann@unibas.ch"
__license__ = "GPL"

"""Create bedGraph for each sample or the whole atlas

This script creates two bedGraph files (for plus and minus strand) from a clusters tsv file. If a sample ID is given, the corresponding tpm column is used; if 'atlas' is given as sample ID, total_tpm is used.
A list of chromosome mappings from ENSEMBL to UCSC names must be specified. 

The input file must contain following columns:
    0 - chr
    1 - strand
    2 - start (1-based coordinate; 1nt will be subtracted here for bed coordinate)
    3 - end
    4 - rep (representative site)
    5 - total_tpm
    6 - annotation
    7 - repSite_signals (";" separated list of polyA signals)
    8 - sample1 tpm
    ...
    X - sampleN tpm

bedGraph output:
    0 - chr
    1 - start
    2 - end
    3 - tpm
"""


# imports
import sys
import os
import time
from argparse import ArgumentParser, RawTextHelpFormatter
import gzip

# Parse input arguments
parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)

parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("-s",
                    "--sample",
                    dest="sample_id",
                    help="sample id or 'atlas'")

parser.add_argument("-i",
                    "--in",
                    dest="input_file",
                    help="input tsv-file containing merged clusters with tpms for all samples (must be gzipped); sample ids have to be present in header")

parser.add_argument("--chr-names",
                    dest="chr_names",
                    nargs="*",
                    help="space separated list of chromosome name pairs that should be included in the output; first part should be the name that occurs in the input file; second part of the pair should be the UCSC equivalent (e.g., 1:chr1)")

parser.add_argument("-p",
                    "--plus",
                    dest="plus_bedgraph",
                    help="output file for plus strand sites")

parser.add_argument("-m",
                    "--minus",
                    dest="minus_bedgraph",
                    help="output file for minus strand sites")

syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):

    # read valid chromosomes for bedgraphs and get the mapping
    valid_chr = {}
    for pair in options.chr_names:
        chrom, ucsc_name = pair.split(":")
        valid_chr[chrom] = ucsc_name

    # If atlas, set sample_id to total_tpm
    if options.sample_id == "atlas":
        sample_id = "total_tpm"
    else:
        sample_id = options.sample_id

    # Handle zipped or unzipped input

    if options.input_file.endswith(".gz"):
        open_input_tsv = gzip.open( options.input_file, "rt")
    else:
        open_input_tsv = open( options.input_file, "r")
    

    ####################################################
    # Loop over all clusters
    # write separate bedgraphs for plus and minus strand
    ####################################################
    
    with open_input_tsv as in_tsv:
        with open(options.plus_bedgraph, "wt") as plus_out, open(options.minus_bedgraph, "wt") as minus_out:

            for line in in_tsv:

                # get the column idx for the sample
                if line.startswith("#"):
                    cols = line.rstrip().split("\t")
                    sample_col = cols.index(sample_id)
                    continue
                
                # get data for each cluster
                F = line.rstrip().split("\t")
                chr    = F[0]
                strand = F[1]
                start  = str(int(F[2])-1)  # subtract 1 for 0-based bed coordinate
                end    = F[3]
                stpm   = str("%.4f" % float(F[sample_col])) # sample tpm
                
                ###################################################
                # Bedgraphs
                # only output lines from valid chromosomes
                # and use ucsc chromosome names
                ###################################################
                bedgraph = []
                # Get the line content for all valid chromosomes
                if chr not in valid_chr:
                    continue
                # Replace ENSEMBL with UCSC chromosome name
                bedgraph.append(valid_chr[chr])
                # We also need start, stop and tpm for sample.
                bedgraph.append(start)
                bedgraph.append(end)
                bedgraph.append(stpm)

                # Depending on strand, write to different output
                if strand == '+':
                    plus_out.write("%s\n" % "\t".join(bedgraph))
                elif strand == '-':
                    minus_out.write("%s\n" % "\t".join(bedgraph))
                else:
                    syserr("No strand information for %s !!\n Ommitting site." %
                           " ".join(bedgraph))


if __name__ == '__main__':
    try:
        try:
            options = parser.parse_args()
        except Exception as e:
            parser.print_help()
            sys.exit()
        if options.verbose:
            start_time = time.time()
            start_date = time.strftime("%d-%m-%Y at %H:%M:%S")
            syserr("############## Started script on %s ##############\n" %
                   start_date)

        main(options)
        if options.verbose:
            syserr("### Successfully finished in %i seconds, on %s ###\n" %
                   (time.time() - start_time,
                    time.strftime("%d-%m-%Y at %H:%M:%S")))
    except KeyboardInterrupt:
        syserr("Interrupted by user after %i seconds!\n" %
               (time.time() - start_time))
        sys.exit(-1)
