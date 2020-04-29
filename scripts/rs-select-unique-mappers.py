__date__ = "2017-02-22"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# imports
import sys
import os
import time
from argparse import ArgumentParser, RawTextHelpFormatter
import gzip

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("-b",
                    "--bed",
                    dest="bed_file",
                    help="input BED-file (must be gzipped)")

syserr = sys.stderr.write
sysout = sys.stdout.write

def main(options):

    if not options.bed_file.endswith(".gz"):
        syserr("[ERROR] input file %s is not gzipped\n"
               % bed_file)
        sys.exit(2)

    read_ids = {}

    # iterate over input file the first time
    with gzip.open(options.bed_file, "rt") as ifile:
        for line in ifile:
            if line.startswith("#"):
                continue
            line_list = line.rstrip().split("\t")
            r_id = line_list[3]
            if r_id in read_ids:
                read_ids[ r_id ] += 1
            else:
                # store read-id, edit distance and counter for this distance
                 read_ids[ r_id ] = 1

    # iterate over the input file once more and only output unique mappers
    with gzip.open(options.bed_file, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                ofile.write(line)
                continue
            l_list = line.rstrip().split("\t")
            r_id = l_list[3]
            if read_ids[ r_id ] == 1:
                # store count (1) as score in the bed column 4
                l_list[3] = "1"
                sysout("%s\n" % "\t".join(l_list))


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
