#!/usr/bin/env python
"""
My template:
"""

__date__ = "2016-07-07"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# imports
import sys
import os
import time
from argparse import ArgumentParser, RawTextHelpFormatter
import shutil
import subprocess
from subprocess import run

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("-s",
                    "--srr_id",
                    dest="srr_id",
                    help="SRR number (from SRA) for the current sample")

parser.add_argument("--outdir",
                    dest="outdir",
                    help="directory to which the fastq files are written")

parser.add_argument("--paired",
                    dest="paired",
                    action="store_true",
                    help="indicate if two samples (paired-end sequencing) belong to that sample id")



# redefine a functions for writing to stdout and stderr to save some writting
syserr = sys.stderr.write
sysout = sys.stdout.write


def main(options):
    """Download the fastq file(s) for the given id"""

    # get the path to the installation of aspera
    starts_with_sep = shutil.which("ascp").startswith(os.sep)
    aspera_path_list = shutil.which("ascp").split(os.sep)
    bin_dir = len(aspera_path_list)
    for i in range(len(aspera_path_list) - 1, -1, -1):
        if aspera_path_list[i] == "bin":
            bin_dir = i
            break

    aspera_path = os.path.join( *aspera_path_list[0:bin_dir] )
    # prepend a directory separator if necessary
    if starts_with_sep and not aspera_path.startswith(os.sep):
        aspera_path = os.sep + aspera_path

    command_list = ["ascp", "-QT", "-l", "300m", "-P33001", "-d", "-i"]
    command_list.append( aspera_path + "/etc/asperaweb_id_dsa.openssh")

    base_address = "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq"

    # SRR, ERR or DRR?
    prefix = options.srr_id[0:3]

    srr_number = options.srr_id.replace(prefix, "")
    if len(srr_number) == 6:
        address = os.path.join(base_address,
                               options.srr_id[:6],
                               options.srr_id,
                               options.srr_id
        )
    elif len(srr_number) == 7:
        address = os.path.join(base_address,
                               options.srr_id[:6],
                               "00" + str(srr_number[-1]),
                               options.srr_id,
                               options.srr_id
        )
    elif len(srr_number) == 8:
        address = os.path.join(base_address,
                               options.srr_id[:6],
                               "0" + str(srr_number[-2:]),
                               options.srr_id,
                               options.srr_id
        )
    elif len(srr_number) == 8:
        address = os.path.join(base_address,
                               options.srr_id[:6],
                               str(srr_number[-3:]),
                               options.srr_id,
                               options.srr_id
        )
    else:
        syserr("[ERROR] SRR id %s has unexpected format. Expected is the form: SRR\d+ with \d+ being 6 to 9 digits\n" % srr_number)
        sys.exit(2)

    if options.paired:
        for read in [1,2]:
            fulladdress = address + "_" + str(read) + ".fastq.gz"
            command = command_list + [fulladdress, options.outdir]
            # print( command )
            returnObj = run(command, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL)
            # syserr("[INFO] command args: %s\n" % str(return_obj.args))
            if returnObj.returncode != 0:
                syserr("[ERROR] command failed\n")
                syserr("[ERROR] command: %s\n" % command)
                sys.exit(2)
    else:
        fulladdress = address + ".fastq.gz"
        command = command_list + [fulladdress, options.outdir]
        return_val = run(command, stderr = subprocess.DEVNULL, stdout = subprocess.DEVNULL).returncode
        if return_val != 0:
            syserr("[ERROR] command failed\n")
            syserr("[ERROR] command: %s\n" % command)
            sys.exit(2)

if __name__ == '__main__':
    try:
        # check if aspera's ascp is available
        if not shutil.which("ascp"):
            syserr("[ERROR] Could not find Aspera's ascp\n")
            syserr("[ERROR] Ensure that ascp is available and rerun the script\n")
            sys.exit(2)
        try:
            options = parser.parse_args()
        except Exception:
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
