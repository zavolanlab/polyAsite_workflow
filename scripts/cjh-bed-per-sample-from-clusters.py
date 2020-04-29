__date__ = "2019-07-24"
__author__ = "Christina J. Herrmann"
__email__ = "christina.herrmann@unibas.ch"
__license__ = "GPL"

"""Create bed for each sample

This script creates a bed file for the specified sample from a clusters tsv file. If sample==atlas, total tpm insted of sample tpm will be written in col 4 (and 8)

The input file must contain following columns:
    0 - chr
    1 - strand
    2 - start (1-based coordinate; 1nt will be subtracted here for bed coordinate)
    3 - end
    4 - rep (representative site)
    5 - total_tpm
    6 - fraction of samples supporting the cluster
    7 - number of protocols supporting the cluster
    8 - annotation
    9 - repSite_signals (";" separated list of polyA signals)
    10 - sample1 tpm
    ...
    X - sampleN tpm

output bed:
    0 - chr
    1 - start
    2 - end 
    3 - id
    4 - sample tpm
    5 - strand
    6 - frac_samples
    7 - nr_prots
    8 - avg tpm accross all samples in atlas
    9 - annotation
    10 - polyA signals
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
                    help="sample id")
parser.add_argument("-i",
                    "--in",
                    dest="input_file",
                    help="input tsv-file containing merged clusters with tpms for all samples; sample ids have to be present in header")
parser.add_argument("-o",
                    "--out",
                    dest="output_bed",
                    help="output bed file with tpm per sample")

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


    # dealing with atlas or sample?
    if options.sample_id == 'atlas':
        sample_id = 'total_tpm'
    else:
        sample_id = options.sample_id

    ####################################################
    # Loop over all clusters
    # write out bedfile for sample
    ####################################################

    with open_input(options.input_file) as in_tsv:
        with open(options.output_bed, "wt") as bed_out:
            
            # Get info from header line
            line = in_tsv.readline()
            # get the column idx for the sample
            try:    
                cols = line.rstrip().split("\t")
                sample_col = cols.index(sample_id)
            except Exception as e:
                syserr("Could not infer sample column.\nCheck whether header line is present.")
                sys.exit()


            # Process remaining lines
            line = in_tsv.readline()
            while line: 
                
                # get data for each cluster
                F = line.rstrip().split("\t")
                chr    = F[0]
                strand = F[1]
                start  = str(int(F[2])-1)  # subtract 1 for 0-based bed coordinate
                end    = F[3]
                rep    = F[4]
                ttpm   = str("%.4f" % float(F[5])) # total tpm
                samp   = F[6]
                prot   = F[7]
                anno   = F[8]
                sign   = F[9]
                stpm   = str("%.4f" % float(F[sample_col])) # sample tpm
                id = ":".join([ chr, rep, strand ])
                
                ####################################################
                # Bed file
                ####################################################
                bed = []
                bed.append(chr) 
                bed.append(start)
                bed.append(end)
                bed.append(id)
                bed.append(stpm) 
                bed.append(strand)
                bed.append(samp)
                bed.append(prot)
                bed.append(ttpm)
                bed.append(anno) 
                bed.append(sign) 
                bed_out.write("%s\n" % "\t".join(bed))

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
