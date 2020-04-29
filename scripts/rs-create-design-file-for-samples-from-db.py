# -*- coding: utf-8 -*-
"""
Created on Tue Jan 03 09:30:36 2017

Add information for TCGA samples to MongoDB

@author: schmiral
"""
__date__ = "Tue Jan 03 09:30:36 2017"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# import
import sys
import os
from argparse import ArgumentParser, RawTextHelpFormatter
import pymongo
import datetime
import time
import json
import yaml
import urllib

# parse input arguments

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")
                    
parser.add_argument("--host",
                    dest="host",
                    help="name of the host for the mongodb")

parser.add_argument("--organism",
                    dest="organism",
                    help="Organism name for which the design file should be created" )

parser.add_argument("--db",
                    dest="db",
                    help="name of the db to query" )

parser.add_argument("--samples-collection",
                    dest="collection",
                    help="name of the collection with the samples" )


syserr = sys.stderr.write
sysout = sys.stdout.write



def time_now():#return time
    curr_time = datetime.datetime.now()
    return curr_time.strftime("%c")

def main(options):

    # parse options
    host = options.host
    organism = options.organism
    collection = options.collection
    db_name = options.db

    ###
    # establish connection to the mongoDB
    ###
    
    # use authentication if user and password are set
    connection = pymongo.MongoClient("mongodb://%s" % host)
    db = connection[ db_name ]
    samples = db[ collection ]

    entries = samples.find({"organism": organism})

    if entries is None:
        syserr("[ERROR] No entries found for %s in %s! Interrupt!\n"
               % (organism, db + ":" + collection) )
        return None

    # iterate over entries and output the required infos
    sysout("sample_id\tSRR\tseries_id\tsource\tsex\tprotocol\treadlen\tfiveAdapter\tthreeAdapter\tpaired\treverse_compl\ttitle\ttreatment\tpubmed\n")

    for entry in entries:
        if not entry["visible"]:
            continue
        sysout("%s\t%s\t%s\t%s\t%s\t%s\t%i\t%s\t%s\t%s\t%r\t%s\t%s\t%s\n"
               % (entry["GEOsampleID"],
                  entry["SRAfile"],
                  entry["GEOseriesID"],
                  entry["source"],
                  entry["sex"],
                  entry["protocol"],
                  entry["reads"]["raw"]["length"]["max"],
                  entry["adaptor5p"],
                  entry["adaptor3p"],
                  entry["pairedEnd"],
                  entry["tobereversed"],
                  entry["title"],
                  entry["treatment"],
                  entry["PubMedID"]
                  )
               )

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

