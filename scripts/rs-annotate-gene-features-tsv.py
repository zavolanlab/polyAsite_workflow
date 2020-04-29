
"""Get annotation features for 3'ends/clusters

    This script is adapted from rs-annotate-gene-features.py, which takes a bed file as input.
    For this script input is .tsv and needs to have at least the following columns:
    0 - chr
    1 - strand
    2 - start; in 1-based coordinate! 1 will be subtracted to pass bed like (0-based) coordinate to function; output will be in 1-based coordinate again.
    3 - end
    4 - x
    5 - x
    
    The annotation will be inserted at index 6. All input columns will be kept (those from index 6 on will be shifted by 1 to the right)
"""

__date__ = "2017-02-14"
__author__ = "Ralf Schmidt"
__email__ = "ralf.schmidt@unibas.ch"
__license__ = "GPL"

# imports
import sys
import os
import time
from argparse import ArgumentParser, RawTextHelpFormatter
import gzip
import re
import bisect


parser = ArgumentParser(description=__doc__, formatter_class=RawTextHelpFormatter)
parser.add_argument("-v",
                    "--verbose",
                    dest="verbose",
                    action="store_true",
                    default=False,
                    help="Be loud!")

parser.add_argument("-g",
                    "--gtf",
                    dest="gtf",
                    help="gene annotation file in uncompressed gtf format")

parser.add_argument("-i",
                    "--input",
                    dest="input_file",
                    help="tsv-file of the poly(A) sites (gzipped or uncompressed; raw 3' ends or clustered poly(A) sites"),

parser.add_argument("--ds-range",
                    dest="ds_range",
                    type = int,
                    help=("number of nt downstream of 3' end and " +
                          "upstream of 5' start (rev. strand) to define the downstream region") )

syserr = sys.stderr.write
sysout = sys.stdout.write

def insert_feature( coordinates_list, sets_list, start, end, term):
    '''
    Insert a feature region(start, end, term) into the
    coordinates_list and sets_list which belong to the features
    objects
    (of the currently processed chromosome and strand)
    '''

    # end is in BED-format -> this means, that it is actually the first
    # pos of the next region

    # set initial values for both lists if necessary
    if len(coordinates_list) == 0:
        coordinates_list.append(0)
        sets_list.append(set())

    inserted = False

    insert_index = bisect.bisect(coordinates_list, start)
    if term not in sets_list[ insert_index - 1 ]:
        if coordinates_list[ insert_index - 1] == start:
            # add term to set of features
            sets_list[ insert_index - 1].add(term)
            # special case because start and coordinate are the same
            insert_index -= 1
            inserted = True
        else:
            # start is bigger than coordinate before
            # create a new step
            coordinates_list.insert(insert_index, start)
            sets_list.insert(insert_index, {term} | sets_list[insert_index - 1] )
            inserted = True

    else:
        # term is already present in previous entry
        # -> no insertion needed
        # but next entry upstream must not be skipped
        insert_index -= 1

    # start-pos handling finished

    # now, iterate over all coordinates that are smaller than the end
    duplicates = set()
    end_index = bisect.bisect( coordinates_list, end )
    for idx in range(insert_index + 1, end_index):
        # special case in last iteration: end and end of current feature are the same
        # in this case: DO NOTHING! because the associated set is not affected by the current feature
        if idx == end_index - 1 and coordinates_list[ idx ] == end:
            break

        if term in sets_list[ idx ]:
            inserted = False
        else:
            sets_list[ idx ].add(term)
            inserted = True
            if sets_list[idx] == sets_list[idx-1]:
                duplicates.add( idx )
            if (idx+1) < len(sets_list):
                if sets_list[idx] == sets_list[ idx+1 ]:
                    duplicates.add(idx+1)

    # handle potential duplicate entries
    for idx in sorted(list(duplicates))[::-1]:
        del coordinates_list[idx]
        del sets_list[idx]

    end_index = bisect.bisect( coordinates_list, end )

    # handle current region end
    # treat special case that end coordinate and current end in list
    # are equal (do nothing in this case)
    if coordinates_list[ end_index -1 ] == end:
        return( coordinates_list, sets_list )

    # insert entry to coordinates_list/sets_list if necessary
    if inserted:
        # term was inserted for the last step
        # insert new entry
        coordinates_list.insert(end_index, end)
        sets_list.insert(end_index, sets_list[end_index - 1].difference({term}) )

    # # test -->
    # syserr("[TMP] Intermediate: coords: %s\n" % str(coordinates_list))
    # syserr("[TMP] Intermediate: sets: %s\n\n" % str(sets_list))
    # # <-- test

    return( coordinates_list, sets_list)

def insert_feature_curr_transcript( coordinates_list, sets_list, start, end, term ):
    '''
    Same function as above but for a single transcript.
    In the case of a single transcript, no leading/trailing empty set is necessary
    '''

    # initialize lists if necessary
    if len(coordinates_list) == 0:
        return( [start,end], [{term},{''}])

    inserted = False

    insert_index = bisect.bisect(coordinates_list, start)

    if term not in sets_list[ insert_index - 1 ]:
        if coordinates_list[ insert_index - 1] == start:
            # add term to set of features
            sets_list[ insert_index - 1].add(term)
            # special case because start and coordinate are the same
            insert_index -= 1
            inserted = True
        else:
            # start is bigger than coordinate before
            # create a new step
            coordinates_list.insert(insert_index, start)
            sets_list.insert(insert_index, {term} | sets_list[insert_index - 1] )
            inserted = True

    else:
        # term is already present in previous entry
        # -> no insertion needed
        # but next entry upstream must not be skipped
        insert_index -= 1

    ### start pos handling finished ###

    # now, iterate over all coordinates that are smaller than the end
    duplicates = set()
    end_index = bisect.bisect( coordinates_list, end )
    for idx in range(insert_index + 1, end_index):
        # special case in last iteration: end and end of current feature are the same
        # in this case: DO NOTHING! because the associated set is not affected by the current feature
        if idx == end_index - 1 and coordinates_list[ idx ] == end:
            break

        if term in sets_list[ idx ]:
            inserted = False
        else:
            sets_list[ idx ].add(term)
            inserted = True
            if sets_list[idx] == sets_list[idx-1]:
                duplicates.add( idx )
            if (idx+1) < len(sets_list):
                if sets_list[idx] == sets_list[ idx+1 ]:
                    duplicates.add(idx+1)

    # handle potential duplicate entries
    for idx in sorted(list(duplicates))[::-1]:
        del coordinates_list[idx]
        del sets_list[idx]

    end_index = bisect.bisect( coordinates_list, end )

    # handle current region end
    # treat special case that end coordinate and current end in list
    # are equal (do nothing in this case)
    if coordinates_list[ end_index -1 ] == end:
        return( coordinates_list, sets_list )

    # insert entry to coordinates_list/sets_list if necessary
    if inserted:
        # term was inserted for the last step
        # insert new entry
        coordinates_list.insert(end_index, end)
        sets_list.insert(end_index, sets_list[end_index - 1].difference({term}) )

    return( coordinates_list, sets_list)


def clean_up_transcript_features(feature_dict, priorities_dict):
    '''
    For a single transcript, only one region annotation is allowed per region.
    For each region, select the feature with the highest priority according
    to the provided dict
    '''

    coordinates_list = feature_dict['coords']
    sets_list = feature_dict['vals']

    # first iteration
    # for each region choose the feature with highest priority
    last_idx = len( sets_list )
    for idx in range(last_idx):
        select_feat = ''
        select_prio = -1
        for feat in sets_list[idx]:
            if feat == '':
                continue
            if priorities_dict[ feat ] > select_prio:
                select_prio = priorities_dict[ feat ]
                select_feat = feat
        # make the selected feature the only one in this region
        sets_list[idx] = select_feat


    # second iteration
    # delete consecutive region when they have the same feature
    duplicates = set()
    for idx in range(1, last_idx):
        if sets_list[idx - 1] == sets_list[idx]:
            duplicates.add(idx)
    for idx in sorted(list(duplicates))[::-1]:
        del coordinates_list[idx]
        del sets_list[idx]

    feature_dict['coords'] = coordinates_list
    feature_dict['vals'] = sets_list

    return feature_dict

def get_all_transcript_features(transcript_features, priorities_dict, ds_range):
    '''
    Iterate over the full dict of features and corresponding regions
    for an entire transcript and assemble the complete lists that contain
    the start coordinate and the feature for each region
    '''

    transcript_feature_dict = {"coords": [], "vals": set()}
    antisense_feature_dict = {"coords": [], "vals": set()}

    # first process the exon entries to get the exon and intron regions
    if "exon" not in transcript_features:
        syserr("[ERROR] No exon found for transcript\n")
        sys.exit(2)

    # order the list of exons to have increasing coordinates
    if (transcript_features['exon'][ 0 ][2] == "-" and
        transcript_features['exon'][ 0 ][3] > transcript_features['exon'][ -1 ][3]):
        # reverse order
        transcript_features['exon'] = transcript_features['exon'][::-1]

    if len(transcript_features['exon']) == 1:
        # process single exon transcripts
        start = transcript_features['exon'][0][3]
        end = transcript_features['exon'][0][4]
        (transcript_feature_dict["coords"],
         transcript_feature_dict["vals"]) = insert_feature_curr_transcript( transcript_feature_dict["coords"],
                                                                            transcript_feature_dict["vals"],
                                                                            start,
                                                                            end,
                                                                            "EX"
         )
    else:
        for ex_idx in range(1,len(transcript_features['exon'])):
            prev_start = transcript_features['exon'][ ex_idx - 1][3]
            prev_end = transcript_features['exon'][ ex_idx - 1][4]
            start = transcript_features['exon'][ ex_idx ][3]
            end = transcript_features['exon'][ ex_idx ][4]
            (transcript_feature_dict["coords"],
             transcript_feature_dict["vals"]) = insert_feature_curr_transcript( transcript_feature_dict["coords"],
                                                                                transcript_feature_dict["vals"],
                                                                                prev_start,
                                                                                prev_end,
                                                                                "EX"
             )
            (antisense_feature_dict["coords"],
             antisense_feature_dict["vals"]) = insert_feature_curr_transcript( antisense_feature_dict["coords"],
                                                                                antisense_feature_dict["vals"],
                                                                                prev_start,
                                                                                prev_end,
                                                                                "AE"
             )
            # store the space inbetween as intron
            (transcript_feature_dict["coords"],
             transcript_feature_dict["vals"]) = insert_feature_curr_transcript( transcript_feature_dict["coords"],
                                                                                transcript_feature_dict["vals"],
                                                                                prev_end,
                                                                                start,
                                                                                "IN"
             )
            (antisense_feature_dict["coords"],
             antisense_feature_dict["vals"]) = insert_feature_curr_transcript( antisense_feature_dict["coords"],
                                                                                antisense_feature_dict["vals"],
                                                                                prev_end,
                                                                                start,
                                                                                "AI"
             )
            if ex_idx == len(transcript_features['exon']) - 1:
                # last exon; process it as well
                (transcript_feature_dict["coords"],
                 transcript_feature_dict["vals"]) = insert_feature_curr_transcript( transcript_feature_dict["coords"],
                                                                                    transcript_feature_dict["vals"],
                                                                                    start,
                                                                                    end,
                                                                                    "AI"
                 )
                (antisense_feature_dict["coords"],
                 antisense_feature_dict["vals"]) = insert_feature_curr_transcript( antisense_feature_dict["coords"],
                                                                                   antisense_feature_dict["vals"],
                                                                                   start,
                                                                                   end,
                                                                                   "AE"
                 )


    # annotate the terminal exon
    # and the downstream regions
    term_ex_idx = -1
    first_ex_idx = 0
    if transcript_features['exon'][ 0 ][2] == "-":
        term_ex_idx = 0
        first_ex_idx = -1

    (transcript_feature_dict["coords"],
     transcript_feature_dict["vals"]) = insert_feature_curr_transcript( transcript_feature_dict["coords"],
                                                                        transcript_feature_dict["vals"],
                                                                        transcript_features['exon'][ term_ex_idx ][3],
                                                                        transcript_features['exon'][ term_ex_idx ][4],
                                                                        "TE"
     )
    if transcript_features['exon'][ 0 ][2] == "+":
        (transcript_feature_dict["coords"],
         transcript_feature_dict["vals"]) = insert_feature_curr_transcript( transcript_feature_dict["coords"],
                                                                            transcript_feature_dict["vals"],
                                                                            transcript_features['exon'][ term_ex_idx ][4],
                                                                            transcript_features['exon'][ term_ex_idx ][4] + ds_range,
                                                                            "DS"
         )
        (antisense_feature_dict["coords"],
         antisense_feature_dict["vals"]) = insert_feature_curr_transcript( antisense_feature_dict["coords"],
                                                                            antisense_feature_dict["vals"],
                                                                            transcript_features['exon'][first_ex_idx][3] - ds_range,
                                                                            transcript_features['exon'][first_ex_idx][3],
                                                                            "AU"
         )
    else:
        # minus strand
        (transcript_feature_dict["coords"],
         transcript_feature_dict["vals"]) = insert_feature_curr_transcript( transcript_feature_dict["coords"],
                                                                            transcript_feature_dict["vals"],
                                                                            transcript_features['exon'][ term_ex_idx ][3] - ds_range,
                                                                            transcript_features['exon'][ term_ex_idx ][3],
                                                                            "DS"
         )
        (antisense_feature_dict["coords"],
         antisense_feature_dict["vals"]) = insert_feature_curr_transcript( antisense_feature_dict["coords"],
                                                                            antisense_feature_dict["vals"],
                                                                            transcript_features['exon'][ first_ex_idx ][4],
                                                                            transcript_features['exon'][ first_ex_idx ][4] + ds_range,
                                                                            "AU"
         )

    # don't allow two features to be annotated at once
    # e.g., delete "EX" where "TE" is annotated
    transcript_feature_dict = clean_up_transcript_features(transcript_feature_dict, priorities_dict)
    antisense_feature_dict = clean_up_transcript_features(antisense_feature_dict, priorities_dict)

    return (transcript_feature_dict, antisense_feature_dict)


def read_gtf(gtf, ds_range, features_with_priorities):
    '''
    iterate through gtf file and return a dict that contains
    the feature for every genomic position
    '''

    feature_dict = {}
    tmp_transcript_dict = {}
    tmp_antisense_dict = {}
    tmp_transcript_id = ''
    tmp_transcript_key = ''
    tmp_antisense_key = ''
    tmp_transcript_features = {}
    tmp_transcript_chrom = ''
    tmp_transcript_strand = ''
    utr_entries = {}

    if gtf.endswith(".gz"):
        in_gtf = gzip.open( gtf, "rt")
    else:
        in_gtf = open( gtf, "r")
    for line in in_gtf:
        if line.startswith("#"):
            continue

        F = line.rstrip().split("\t")

        if F[2] == "transcript":
            # get transcript id
            mo = re.match('.+transcript_id\s\"([^\"]+)', F[8])
            assert mo
            tr_id = mo.groups()[0]
            chrom = F[0]
            strand = F[6]
            # use chr:strand as key
            key = chrom + ":" + strand
            antisense_key = chrom + ":+"
            if strand == "+":
                antisense_key = chrom + ":-"
            if key not in feature_dict:
                feature_dict[ key ] = {"coords": [], "sets": []}
            if antisense_key not in feature_dict:
                feature_dict[ antisense_key ] = {"coords": [], "sets": []}

            if tmp_transcript_id != '':
                # process previous transcript
                (tmp_transcript_dict,
                 tmp_antisense_dict)= get_all_transcript_features(tmp_transcript_features,
                                                                  features_with_priorities,
                                                                  ds_range
                 )

                for idx in range(len(tmp_transcript_dict["coords"])-1):
                    (feature_dict[ tmp_transcript_key ]["coords"],
                     feature_dict[ tmp_transcript_key ]["sets"]) = insert_feature( feature_dict[ tmp_transcript_key ]["coords"],
                                                                                 feature_dict[ tmp_transcript_key ]["sets"],
                                                                                 tmp_transcript_dict["coords"][idx],
                                                                                 tmp_transcript_dict["coords"][idx+1],
                                                                                 tmp_transcript_dict["vals"][idx]
                     )
                for idx in range(len(tmp_antisense_dict["coords"])-1):
                    (feature_dict[ tmp_antisense_key ]["coords"],
                     feature_dict[ tmp_antisense_key]["sets"]) = insert_feature( feature_dict[ tmp_antisense_key ]["coords"],
                                                                                 feature_dict[ tmp_antisense_key ]["sets"],
                                                                                 tmp_antisense_dict["coords"][idx],
                                                                                 tmp_antisense_dict["coords"][idx+1],
                                                                                 tmp_antisense_dict["vals"][idx]
                     )

            tmp_transcript_dict = {}
            tmp_antisense_dict = {}
            tmp_transcript_features = {}
            tmp_transcript_chrom = chrom
            tmp_transcript_strand = strand
            tmp_transcript_id = tr_id
            tmp_transcript_key = key
            tmp_antisense_key = antisense_key


        #######################################################
        #######################################################

        elif F[2] == "exon":
            # get transcript id
            mo = re.match('.+transcript_id\s\"([^\"]+)', F[8])
            assert mo
            tr_id = mo.groups()[0]
            chrom = F[0]
            strand = F[6]
            # change coordinates to BED format
            start = int(F[3]) - 1
            end = int(F[4])

            # consistency check
            if tr_id != tmp_transcript_id:
                syserr("[ERROR] Transcript id of current exon does not match the currently processed transcript\n")
                syserr("[ERROR] Exon transcript: %s, current transcript: %s\n"
                       % (tr_id, tmp_transcript_id))
                sys.exit(2)

            if 'exon' not in tmp_transcript_features:
                tmp_transcript_features['exon'] = []
            exon_infos = [tr_id, chrom, strand, start, end]
            tmp_transcript_features['exon'].append( exon_infos )

    in_gtf.close()
    ###
    # process last transcript
    ###
    if tmp_transcript_id != '':
        # process previous transcript
        (tmp_transcript_dict,
         tmp_antisense_dict)= get_all_transcript_features(tmp_transcript_features,
                                                          features_with_priorities,
                                                          ds_range
        )

        for idx in range(len(tmp_transcript_dict["coords"])-1):
            (feature_dict[ tmp_transcript_key ]["coords"],
             feature_dict[ tmp_transcript_key ]["sets"]) = insert_feature( feature_dict[ tmp_transcript_key ]["coords"],
                                                                           feature_dict[ tmp_transcript_key ]["sets"],
                                                                           tmp_transcript_dict["coords"][idx],
                                                                           tmp_transcript_dict["coords"][idx+1],
                                                                           tmp_transcript_dict["vals"][idx]
             )

    return(feature_dict)

def get_overlapping_features( feat_dict, start, end):
    '''
    returns a set of features that overlaps
    with the segment defined by start and end
    '''

    coords_list = feat_dict[ 'coords' ]
    feature_sets_list = feat_dict[ 'sets' ]

    return_set = set()

    start_idx = bisect.bisect( coords_list, start) - 1
    end_idx = bisect.bisect( coords_list, end)

    if start_idx == end_idx:
        end_idx += 1

    for region_idx in range( start_idx, end_idx ):
        if region_idx == end_idx - 1:
            # only count this features if
            # the corresponding coordinate is << then
            # the end coordinate (BED format: end coord is not included in the feature)
            if end ==  coords_list[ region_idx ]:
                continue
        return_set = return_set.union( feature_sets_list[ region_idx ] )

    return( return_set )


def main(options):
    '''
    read the gtf_file
    and iterate over the input features to
    define their overlap with the features given in the gtf
    '''

    # The following features are used:
    features_with_priorities = {}
    prio = 11
    for f in ["TE", "EX", "IN", "DS", "AE", "AI", "AU"]:
        features_with_priorities[ f ] = prio
        prio -= 1

    feature_dict = read_gtf( options.gtf, options.ds_range, features_with_priorities )

    if options.verbose:
        syserr("[INFO] Finished reading gtf. Start iterating over the input file\n")

    ###########################
    # read input tsv file
    # and insert the annotation as new 7th column
    ###########################

    if options.input_file.endswith(".gz"):
        open_input_tsv = gzip.open( options.input_file, "rt")
    else:
        open_input_tsv = open( options.input_file, "r")

    with open_input_tsv as input_tsv:
        for line in input_tsv:
    
            # print header
            if line.startswith("#"):
                line_list = line.rstrip().split("\t")
                line_list.insert(6,"annotation")
                sysout("%s\n" % "\t".join(line_list))
                continue
    
    
            line_list = line.rstrip().split("\t")
            curr_key = line_list[0] + ":" + line_list[1]
            curr_start = int(line_list[2])-1 # subtract 1nt from start for bed coordinate
            curr_end = int(line_list[3])
    
            if curr_key not in feature_dict:
                # insert the annotation as new 7th column
                line_list.insert(6,"IG")
                sysout("%s\n" % "\t".join(line_list))
                continue
    
            feature_set = get_overlapping_features( feature_dict[ curr_key ], curr_start, curr_end)
    
            number_of_features = len(feature_set)
            if (number_of_features == 0 or
                ( number_of_features == 1 and "" in feature_set)):
                # insert the annotation as new 7th column
                line_list.insert(6,"IG")
                sysout("%s\n" % "\t".join(line_list))
            else:
                max_prio = 0
                selected_feature = None
                feature_set.discard("")
                for feat in feature_set:
                    if features_with_priorities[ feat ] > max_prio:
                        max_prio = features_with_priorities[ feat ]
                        selected_feature = feat
                # insert the annotation as new 7th column
                line_list.insert(6,selected_feature)
                sysout("%s\n" % "\t".join(line_list))
    

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
