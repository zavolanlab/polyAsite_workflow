__date__ = "2019-08-06"
__author__ = "Christina J. Herrmann"
__email__ = "christina.herrmann@unibas.ch"
__license__ = "GPL"

"""Calculate sample support for clusters

This script calculates for each cluster in input clusters.anno.tsv.gz the fraction
of samples, as well as the number of protocols that support the cluster.

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

output tsv:
    0 - chr
    1 - strand
    2 - start (1-based coordinate; 1nt will be subtracted here for bed coordinate)
    3 - end
    4 - rep (representative site)
    5 - total_tpm
    6 - fraction of supporting samples
    7 - number of supporting protocols
    8 - annotation
    9 - repSite_signals (";" separated list of polyA signals)
    10 - sample1 tpm
    ...
    X - sampleN tpms
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
parser.add_argument("-d",
                    "--design",
                    dest="design_file",
                    help="design file (tsv)")
parser.add_argument("-i",
                    "--in",
                    dest="input_file",
                    help="input tsv-file containing merged clusters with tpms for all samples; sample ids have to be present in header")
parser.add_argument("-o",
                    "--out",
                    dest="output_file",
                    help="output tsv file.")

syserr = sys.stderr.write
sysout = sys.stdout.write


def open_input(filename):
    """Handle zipped or unzipped input"""
    if filename.endswith(".gz"):
        handler = gzip.open( filename, "rt")
    else:
        handler = open( filename, "r")
    return handler


def main(options):

    protocols = []
    samples = []
    ####################################################
    # Get protocol for each sample
    ####################################################
    with open(options.design_file, "r") as design:
        # Read the header line to get index of protocol column
        line = design.readline()
        cols = line.rstrip().split("\t")
        id_col = cols.index("sample_id")
        prot_col = cols.index("protocol")

        # Proceed to the remaining lines
        line = design.readline().rstrip().split("\t")
        while len(line) > 1:
            samples.append(line[id_col])
            protocols.append(line[prot_col])
            line = design.readline().rstrip().split("\t")

    syserr("Finished reading design file.\n")
    syserr(" ".join(samples))
    syserr(" ".join(protocols))
    ####################################################
    # Loop over all clusters
    ####################################################

    with open_input(options.input_file) as in_tsv:
        with open(options.output_file, "wt") as out:

            #############################
            # Header
            ############################ 
            # Get index info from header
            line = in_tsv.readline()
            cols = line.rstrip().split("\t")
            try:
                first_sample_col = cols.index(samples[0])
            except Exception as e:
                syserr("Could not infer first sample column.\nCheck whether header line is present and design file matches cluster file.")
                sys.exit()
            
            # Write new header
            header = cols[0:6]
            header.append("frac_samples")
            header.append("nr_prots")
            header += cols[6:]
            out.write("%s\n" % "\t".join(header))
 
            syserr("Processed header line, proceeding to cluster lines.\n")    
            syserr("First sample column: %s\n" % first_sample_col)

            line = in_tsv.readline()
            while line:
                # get base data for each cluster
                cluster_line = line.rstrip().split("\t")

                # calculate support
                # supporting samples
                supsam = 0
                # supporting protocols
                suppro = []

                for s in range(first_sample_col, len(cluster_line)):
                    if float(cluster_line[s]) > 0:
                        supsam += 1
                        suppro.append(protocols[s-first_sample_col])
                syserr("supporting samples %i\n" % int(supsam))
                supsams = "{0:.2f}".format( supsam / len(samples))
                syserr("supporting samples fraction %s\n" % str(supsams))
                supprots = str(len(set(suppro)))



                ####################################################
                # Output file
                ####################################################
                cluster_out = []
                cluster_out += cluster_line[0:6]
                cluster_out.append(supsams)
                cluster_out.append(supprots)
                cluster_out += cluster_line[6:]

                out.write("%s\n" % "\t".join(cluster_out))
                
                line = in_tsv.readline()



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
