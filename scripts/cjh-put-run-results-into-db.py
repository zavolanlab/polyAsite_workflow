# -*- coding: utf-8 -*-
"""
Created on Thu May 31 09:30:36 2018

Add information of a run to prepare
a poly(A) site atlas to the MongoDB

@author: schmiral
"""
__date__ = "Thu May 31 09:30:36 2018"
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
import pandas as pd

# parse input arguments

parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("--config",
                    dest="config",
                    help=("file name providing the config " +
                        "information of the run (yaml format required)") )

parser.add_argument("--design-file",
                    dest="design_file",
                    help=("tsv file with sample specific information; one sample per row with ids in first column") )


syserr = sys.stderr.write
sysout = sys.stdout.write

def obtain_infos(sampleID, samples, stats, config_dict):
    '''
    include the information for a new sample
    '''

    data_dict = {}

    # include all stats stored during the run
    for entry in stats:
        F = entry.split(".")
        if "raw" in F:
            # the entry can be put in the db as is
            if len(F) == 3:
                if F[0] not in data_dict:
                    data_dict[ F[0] ] = {}
                if F[1] not in data_dict[ F[0] ]:
                    data_dict[F [0] ][F[1] ] = {}
                data_dict[ F[0] ][ F[1] ][ F[2] ] = stats[entry]
            elif len(F) == 4:
                if F[0] not in data_dict:
                    data_dict[ F[0] ] = {}
                if F[1] not in data_dict[ F[0] ]:
                    data_dict[F [0] ][F[1] ] = {}
                if F[2] not in data_dict[ F[0] ][F[1]]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ] = stats[entry]
            elif len(F) == 5:
                if F[0] not in data_dict:
                    data_dict[ F[0] ] = {}
                if F[1] not in data_dict[ F[0] ]:
                    data_dict[F [0] ][F[1] ] = {}
                if F[2] not in data_dict[ F[0] ][F[1]]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                if F[3] not in data_dict[ F[0] ][F[1]][ F[2] ]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ] = {}
                data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ F[4] ] = stats[entry]
        else:
            # the entry contains information that depend on the processing
            # and on the genome -> include this information additionally
            genome = config_dict['genome']
            if len(F) == 3:
                if F[0] not in data_dict:
                    data_dict[ F[0] ] = {}
                if F[1] not in data_dict[ F[0] ]:
                    data_dict[ F[0] ][ F[1] ] = {}
                    data_dict[ F[0] ][ F[1] ][ genome ] = {}
                data_dict[ F[0] ][ F[1] ][ genome ][ F[2] ] = stats[entry]
            elif len(F) == 4:
                if F[0] not in data_dict:
                    data_dict[ F[0] ] = {}
                if F[1] not in data_dict[ F[0] ]:
                    data_dict[F [0] ][F[1] ] = {}
                if F[2] not in data_dict[ F[0] ][F[1]]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                    data_dict[ F[0] ][ F[1] ][ F[2] ][ genome]  = {}
                data_dict[ F[0] ][ F[1] ][ F[2] ][ genome][ F[3] ] = stats[entry]
            elif len(F) == 5:
                if F[0] not in data_dict:
                    data_dict[ F[0] ] = {}
                if F[1] not in data_dict[ F[0] ]:
                    data_dict[F [0] ][F[1] ] = {}
                if F[2] not in data_dict[ F[0] ][F[1]]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                if F[3] not in data_dict[ F[0] ][F[1]][ F[2] ]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ] = {}
                if genome not in data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ genome ] = {}
                data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ genome ][ F[4] ] = stats[entry]
            elif len(F) == 6:
                if F[0] not in data_dict:
                    data_dict[ F[0] ] = {}
                if F[1] not in data_dict[ F[0] ]:
                    data_dict[F [0] ][F[1] ] = {}
                if F[2] not in data_dict[ F[0] ][F[1]]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                if F[3] not in data_dict[ F[0] ][F[1]][ F[2] ]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ] = {}
                if F[4] not in data_dict[ F[0] ][F[1]][ F[2] ][ F[3] ]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ F[4] ] = {}
                if genome not in data_dict[ F[0] ][F[1]][ F[2] ][ F[3] ][ F[4] ]:
                    data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ F[4] ][ genome ] = {}
                data_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ F[4] ][ genome ][ F[5] ] = stats[entry]

    # include additional information
    data_dict["GEOsampleID"] = sampleID
    data_dict["GEOseriesID"] = samples.loc[sampleID, "series_id"]
    data_dict["organism"] = config_dict['organism_name_db']
    data_dict["source"] = samples.loc[sampleID, "source"]
    data_dict["protocol"] = samples.loc[sampleID, "protocol"]
    data_dict["sex"] = samples.loc[sampleID, "sex"]
    data_dict["visible"] = True
    data_dict["treatment"] = samples.loc[sampleID, "treatment"]
    data_dict["title"] = samples.loc[sampleID, "title"]
    data_dict["PubMedID"] = str(samples.loc[sampleID, "pubmed"])

    # get the restrictions for reads to be included
    if "valid" not in data_dict["reads"]:
        data_dict["reads"]["valid"] = {"cutoffs": { } }
    if "cutoffs" not in data_dict["reads"]["valid"]:
        data_dict["reads"]["valid"]["cutoffs"] = {genome: {}}
    data_dict["reads"]["valid"]["cutoffs"][ genome ]["maxAs"] = config_dict['maxAcontent']
    data_dict["reads"]["valid"]["cutoffs"][ genome ]["maxNs"] = config_dict['maxN']
    if samples.loc[sampleID, "protocol"] == "A-seq2":
        data_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - 7
    elif samples.loc[sampleID, "protocol"] == "3'READS":
        data_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - 6
    elif (samples.loc[sampleID, "protocol"] == "A-seq" or
          samples.loc[sampleID, "protocol"] == "3'-Seq (Mayr)"):
        data_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - int(config_dict['min_sense_strand_shortening'])
    elif samples.loc[sampleID, "protocol"] == "3P-Seq":
        data_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - 2
    else:
        data_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"])

    # get information to allow visualization of the 3'ends in the UCSC genome browser
    if "ucsc" not in data_dict:
        data_dict["ucsc"] = {}
    if genome not in data_dict["ucsc"]:
        data_dict["ucsc"][genome] = {}
    data_dict["ucsc"][genome]["organism"] = config_dict["ucsc_organism"]
    data_dict["ucsc"][genome]["db"] = config_dict["ucsc_db"]

    return data_dict

def update_information(sample_dict, samples, stats, config_dict):
    '''
    include new/additional information from the latest run
    or overwrite existing information for an already existing
    sample
    '''
    genome = config_dict['genome']
    sampleID = sample_dict["GEOsampleID"]
    for entry in stats:
        F = entry.split(".")
        if "raw" in F:
            # the entry should be there already
            # skip it
            continue
        else:
            if len(F) == 3:
                if F[0] not in sample_dict:
                    sample_dict[ F[0] ] = {}
                if F[1] not in sample_dict[ F[0] ]:
                    sample_dict[F [0] ][F[1] ] = {}
                if genome not in sample_dict[F [0] ][F[1] ]:
                    sample_dict[F [0] ][F[1] ][ genome ] = {}
                sample_dict[ F[0] ][ F[1] ][ genome ][ F[2] ] = stats[entry]
            elif len(F) == 4:
                if F[0] not in sample_dict:
                    sample_dict[ F[0] ] = {}
                if F[1] not in sample_dict[ F[0] ]:
                    sample_dict[F [0] ][F[1] ] = {}
                if F[2] not in sample_dict[ F[0] ][F[1]]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                if genome not in sample_dict[ F[0] ][ F[1] ][ F[2] ]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ][ genome]  = {}
                sample_dict[ F[0] ][ F[1] ][ F[2] ][ genome][ F[3] ] = stats[entry]
            elif len(F) == 5:
                if F[0] not in sample_dict:
                    sample_dict[ F[0] ] = {}
                if F[1] not in sample_dict[ F[0] ]:
                    sample_dict[F [0] ][F[1] ] = {}
                if F[2] not in sample_dict[ F[0] ][F[1]]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                if F[3] not in sample_dict[ F[0] ][F[1]][ F[2] ]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ] = {}
                if genome not in sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ genome ] = {}
                sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ genome ][ F[4] ] = stats[entry]
            elif len(F) == 6:
                if F[0] not in sample_dict:
                    sample_dict[ F[0] ] = {}
                if F[1] not in sample_dict[ F[0] ]:
                    sample_dict[F [0] ][F[1] ] = {}
                if F[2] not in sample_dict[ F[0] ][F[1]]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ] = {}
                if F[3] not in sample_dict[ F[0] ][F[1]][ F[2] ]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ] = {}
                if F[4] not in sample_dict[ F[0] ][F[1]][ F[2] ][ F[3] ]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ F[4] ] = {}
                if genome not in sample_dict[ F[0] ][F[1]][ F[2] ][ F[3] ][ F[4] ]:
                    sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ F[4] ][ genome ] = {}
                sample_dict[ F[0] ][ F[1] ][ F[2] ][ F[3] ][ F[4] ][ genome ][ F[5] ] = stats[entry]

    # get the restrictions for reads to be included
    if "valid" not in sample_dict["reads"]:
        sample_dict["reads"]["valid"] = {"cutoffs": { } }
    if "cutoffs" not in sample_dict["reads"]["valid"]:
        sample_dict["reads"]["valid"]["cutoffs"] = {genome: {}}
    if genome not in sample_dict["reads"]["valid"]["cutoffs"]:
        sample_dict["reads"]["valid"]["cutoffs"][genome] = {}
    sample_dict["reads"]["valid"]["cutoffs"][ genome ]["maxAs"] = config_dict['maxAcontent']
    sample_dict["reads"]["valid"]["cutoffs"][ genome ]["maxNs"] = config_dict['maxN']
    if samples.loc[sampleID, "protocol"] == "A-seq2":
        sample_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - 7
    elif samples.loc[sampleID, "protocol"] == "3'READS":
        sample_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - 6
    elif (samples.loc[sampleID, "protocol"] == "A-seq" or
          samples.loc[sampleID, "protocol"] == "3'-Seq (Mayr)"):
        sample_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - int(config_dict['min_sense_strand_shortening'])
    elif samples.loc[sampleID, "protocol"] == "3P-Seq":
        sample_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"]) - 2
    else:
        sample_dict["reads"]["valid"]["cutoffs"][ genome ]["maxLength"] = int(samples.loc[sampleID, "readlen"])

    # get information to allow visualization of the 3'ends in the UCSC genome browser
    if "ucsc" not in sample_dict:
        sample_dict["ucsc"] = {}
    if genome not in sample_dict["ucsc"]:
        sample_dict["ucsc"][genome] = {}
    sample_dict["ucsc"][genome]["organism"] = config_dict["ucsc_organism"]
    sample_dict["ucsc"][genome]["db"] = config_dict["ucsc_db"]

    # no return is necessary because the sample_entry was
    # edited in place and is stored as such again

def obtain_clusters_infos(samples, stats, config_dict):
    '''
    insert a new atlas into the db
    '''

    genome = config_dict["genome"]

    data_dict = {}
    ### --> CJH20190612: Get the stats
    for entry in stats:
        F = entry.split(".")
        if len(F) == 1:
            if F[0] not in data_dict:
                data_dict[ F[0] ] = {}
                data_dict[ F[0] ][genome] = stats[entry]
        elif len(F) == 2:
            if F[0] not in data_dict:
                data_dict[ F[0] ] = {}
            if F[1] not in data_dict[ F[0] ]:
                data_dict[ F[0] ][ F[1] ] = {}
                data_dict[ F[0] ][ F[1] ][genome] = stats[entry]
        elif len(F) == 3:
            if F[0] not in data_dict:
                data_dict[ F[0] ] = {}
            if F[1] not in data_dict[ F[0] ]:
                data_dict[F [0] ][F[1] ] = {}
                data_dict[F [0] ][F[1] ][ genome ] = {}
            data_dict[ F[0] ][ F[1] ][ genome ][ F[2] ] = stats[entry]
    ### <--- CJH20190612

    data_dict["organism"] = config_dict["organism_name_db"]
    data_dict["publicName"] = config_dict["atlas.release_name"]
    data_dict["version"] = config_dict["atlas.release_name"]
    data_dict["releaseDate"] = datetime.datetime.now()
    data_dict['public'] = True
    data_dict["GEOsampleIDs"] = samples.index.tolist()
    data_dict["polyAsignals"] = config_dict["polyA_signals"]
    data_dict["genomeVersion"] = config_dict["genome"]
    # ucsc relevant info to allow display of tracks
    data_dict["ucsc"] = {}
    if genome not in data_dict["ucsc"]:
        data_dict["ucsc"][genome] = {}
    data_dict["ucsc"][genome]["organism"] = config_dict["ucsc_organism"]
    data_dict["ucsc"][genome]["db"] = config_dict["ucsc_db"]

    return data_dict

def update_clusters_information(cluster_dict, samples, stats, config_dict):
    '''
    update information for an existing poly(A) atlas
    '''

    genome = config_dict["genome"]
    for entry in stats:
        F = entry.split(".")
        if len(F) == 1:
            if F[0] not in cluster_dict:
                cluster_dict[ F[0] ] = {}
                cluster_dict[ F[0] ][genome] = stats[entry]
        elif len(F) == 2:
            if F[0] not in cluster_dict:
                cluster_dict[ F[0] ] = {}
            if F[1] not in cluster_dict[ F[0] ]:
                cluster_dict[ F[0] ][ F[1] ] = {}
                cluster_dict[ F[0] ][ F[1] ][genome] = stats[entry]
        elif len(F) == 3:
            if F[0] not in cluster_dict:
                cluster_dict[ F[0] ] = {}
            if F[1] not in cluster_dict[ F[0] ]:
                cluster_dict[F [0] ][F[1] ] = {}
                cluster_dict[F [0] ][F[1] ][ genome ] = {}
            cluster_dict[ F[0] ][ F[1] ][ genome ][ F[2] ] = stats[entry]
    ### <--- CJH20190612


    cluster_dict["releaseDate"] = datetime.datetime.now()
    cluster_dict['public'] = True
    cluster_dict["GEOsampleIDs"] = samples.index.tolist()
    cluster_dict["polyAsignals"] = config_dict["polyA_signals"]
    cluster_dict["genomeVersion"] = config_dict["genome"]

    if "ucsc" not in cluster_dict:
        cluster_dict["ucsc"] = {}
    if genome not in cluster_dict["ucsc"]:
        cluster_dict["ucsc"][genome] = {}
    cluster_dict["ucsc"][genome]["organism"] = config_dict["ucsc_organism"]
    cluster_dict["ucsc"][genome]["db"] = config_dict["ucsc_db"]

    # no return necessary

def update_samples_in_db(samples, config_dict):
    '''
    include all relevant information about
    individual samples in the db
    '''

    # establish connection to the mongoDB

    # use authentication if user and password are set
    if( config_dict['mongoDB.user'] is not None and
        config_dict['mongoDB.password'] is not None):
        pswd = urllib.quote_plus(config_dict['password'])
        connection = pymongo.MongoClient("mongodb://%s:%s@%s"
                                         % (config_dict['user'], pswd, config_dict['mongoDB.host'] ) )
    else:
        connection = pymongo.MongoClient("mongodb://%s"
                                         % config_dict['mongoDB.host'])
    db = connection[ config_dict['mongoDB.db'] ]
    samples_collection = db[ config_dict['mongoDB.samples_collection'] ]

    # iterate over samples of the current run
    for sample_id in samples.index:
        # get the run related information
        numbers_file = os.path.join(config_dict['samples_dir'],
                                    "counts",
                                    (sample_id +
                                     "_" +
                                     config_dict["genome"] +
                                     ".noBG3pSites.out"))
        # consistency check
        if not os.path.isfile(numbers_file):
            syserr("[ERROR] No numbers found for %s under %s; maybe the run is corrupted?\n"
                   % (sample_id, numbers_file) )
            sys.exit(2)

        sample_run_stats = {}
        with open(numbers_file, "r") as numbers_in:
            for line in numbers_in:
                F = line.rstrip().split("\t")
                sample_run_stats[ F[0] ] = int(F[1])

        # check if sample is in db already
        if samples_collection.find_one({"GEOsampleID": sample_id}) != None:
            syserr("[INFO] Already found entry for sample %s. Sample information is updated!\n"
                   % sample_id)

            entry_cursor = samples_collection.find({"GEOsampleID": sample_id})
            if entry_cursor.count() != 1:
                syserr("[ERROR] No unambiguous entry found for %s! Interrupt!\n"
                       % sample_id)
                return None

            entry = list(entry_cursor)[0]
            # get updated information
            update_information(entry,
                               samples,
                               sample_run_stats,
                               config_dict
            )
            samples_collection.save(entry)

        else:
            syserr("[INFO] Insert entry for sample %s.\n"
                   % sample_id)
            data_dict = obtain_infos(sample_id,
                                     samples,
                                     sample_run_stats,
                                     config_dict
            )
            samples_collection.insert_one(data_dict)

def update_clusters_in_db(samples, config_dict):
    '''
    include the information from the created atlas into the db
    '''

    # establish connection to the mongoDB

    # use authentication if user and password are set
    if( config_dict['mongoDB.user'] is not None and
        config_dict['mongoDB.password'] is not None):
        pswd = urllib.quote_plus(config_dict['password'])
        connection = pymongo.MongoClient("mongodb://%s:%s@%s"
                                         % (config_dict['user'], pswd, config_dict['mongoDB.host'] ) )
    else:
        connection = pymongo.MongoClient("mongodb://%s"
                                         % config_dict['mongoDB.host'])
    db = connection[ config_dict['mongoDB.db'] ]
    clusters_collection = db[ config_dict['mongoDB.atlas_collection'] ]

    atlas_organism = config_dict["organism_name_db"]
    genome = config_dict["genome"]
    release_name = config_dict["atlas.release_name"]

    ### --> CJH 20190612: Let's also get some cluster stats
    numbers_file = os.path.join( config_dict['atlas_dir'],
                       config_dict['organism'],
                       config_dict['genome'],
                       config_dict['atlas.release_name'],
                                     "counts",
                                     "clusters.stats.out" )
    # consistency check
    if not os.path.isfile(numbers_file):
        syserr("[ERROR] No numbers found for %s under %s; maybe the run is corrupted?\n"
               % (config_dict["atlas.release_name"], numbers_file) )
        sys.exit(2)

    cluster_run_stats = {}
    with open(numbers_file, "r") as numbers_in:
        for line in numbers_in:
            F = line.rstrip().split("\t")
            cluster_run_stats[ F[0] ] = int(F[1])
    ### <-- CJH20190612

    # check if sample is in db already
    if clusters_collection.find_one({"organism": atlas_organism, "publicName": release_name}) != None:
        syserr("[INFO] Already found entry for atlas %s in %s. Information is updated!\n"
               % (release_name, atlas_organism) )

        entry_cursor = clusters_collection.find({"organism": atlas_organism, "publicName": release_name})
        if entry_cursor.count() != 1:
            syserr("[ERROR] No unambiguous entry found for %s in %s! Interrupt!\n"
                   % (release_name, atlas_organism) )
            return None

        entry = list(entry_cursor)[0]
        # get updated information
        update_clusters_information(entry,  samples, cluster_run_stats, config_dict)
        clusters_collection.save(entry)
    else:
        syserr("[INFO] Insert entry for %s in %s.\n"
               % (release_name, atlas_organism) )
        data_dict = obtain_clusters_infos(samples, cluster_run_stats, config_dict)
        clusters_collection.insert_one(data_dict)


def time_now():
    '''
    return current time
    '''

    curr_time = datetime.datetime.now()
    return curr_time.strftime("%c")

def main(options):

    #----------
    # read the samples information
    #----------
    samples = pd.read_table(options.design_file, index_col = 0, comment = "#")

    #----------
    # parse the yaml file with all config information
    #----------
    with open(options.config) as f:
        config_dict = yaml.load(f)

    ###########

    #----------
    # update/create the individual entries for the samples
    #----------

    update_samples_in_db(samples, config_dict)

    #----------
    # update/create the entry for the poly(A) sites
    #----------

    update_clusters_in_db(samples, config_dict)


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
