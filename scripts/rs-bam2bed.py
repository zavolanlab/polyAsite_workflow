import sys
import os
import optparse
import re
import itertools
import multiprocessing
import time

max_processors = 1
#def __init__(self):

parser = optparse.OptionParser("usage: %prog -p [processors, optional]")
parser.add_option("-p", "--processors", dest="processors", type="int", help="Number of processors to use", default=1)

(options, args) = parser.parse_args(sys.argv)
if None in options.__dict__.values():
    parser.error("Incorrect number of arguments. Did you specify the number of cores?")
    sys.exit(-1)
else:
    if options.processors >= 1:
        max_processors = options.processors


def align():
    ts = time.time()
    ### Multiprocessing ### 
    pool = multiprocessing.Pool(max_processors)
    ### THIS WAS MULTIPROCESSING ###

    #sam_tuples = readfile(input_sam, "\t")
    count = 0
    data_lines = []
    batch_size = 100000

    for sam_line in sys.stdin:
        if sam_line[0] == "@":
            continue
        sam_tuple = tuple(sam_line.split("\t"))
        data_lines.append(sam_tuple)
        if count == batch_size:
            res = pool.map(generate_bed_line, data_lines, max_processors)
            for r in res:
                print(r)
            res = None
            count = 0
            data_lines = []
        else:
            count += 1
    if count > 0:
        res = pool.map(generate_bed_line, data_lines, max_processors)
        for r in res:
            print(r)

def generate_bed_line(sam):

    if sam[:1] == "@":
        return

    for i in range(len(sam)):
        if sam[i].startswith("NM:"):
            # NM stores the edit distance between read and ref
            sam_index_space_pos = i
        elif sam[i].startswith("MD:"):
            # MD stores the string for mismatching positions
            sam_mismatches_pos = i
        elif sam[i].startswith("NH:"):
            # store number of mappings for that read
            read_weight = 1.0 / float(sam[i].replace("NH:i:", ""))
    
    
    # CONSTANTS
    BED_CHR=0
    BED_START=1
    BED_STOP=2
    BED_READ=3
    BED_SCORE=4
    BED_STRAND=5
    BED_MID=6

    SAM_READ_ID=0
    SAM_FLAG=1
    SAM_CHR=2
    SAM_START=3
    SAM_SCORE = 4
    SAM_MID = 5 # aka cigar string
    SAM_STAR = 6 # next entry on same chromosome?
    SAM__ = 7 # mate info
    SAM__ = 8 # mate info
    SAM_READ=9
    SAM_STAR = 10
    SAM__ = 42
    SAM_INDEX_SPACE = sam_index_space_pos
    SAM_MISMATCHES = sam_mismatches_pos

    # 5th digit in SAM_FLAG in binary format tells the strand: 0 -> plus, 1 -> minus

    #Now reading in the sam_tuples sequentially.
    #Need to convert this in some kind of bed format

    MID_length = get_match_length(sam[SAM_MID])

    strand = "?"

    SAM_FLAG_BINARY = '{0:012b}'.format(int(sam[SAM_FLAG]))
    SAM_STRAND = int(SAM_FLAG_BINARY[::-1][4])
    if SAM_STRAND == 0:
        strand = "+"
    elif SAM_STRAND == 1:
        strand = "-"
    else:
        sys.write("[ERROR] Could not infer strand from " + SAM_FLAG + "\n")
        sys.exit(2)

    bed_line = ""
    bed_line += sam[SAM_CHR]+"\t"
    bed_line += str(int(sam[SAM_START])-1)+"\t"
    bed_line += str(int(sam[SAM_START])-1+MID_length)+"\t"
    bed_line += sam[SAM_READ_ID]+"\t"
    # the read weight is later on overwritten by the edit distance
    bed_line += str(read_weight) +"\t"
    bed_line += strand


    #At this point, the bed file is almost complete, just needs the MMID
    #Make a tuple with it!

    bed = bed_line.split("\t")             


    g_length = MID_length

    # genome is set as only "?" because we don't know the seqeuence
    g = "".ljust(g_length, "?")
    g = list(g)

    # get the read that was aligned
    soft_clipped_start = 0
    soft_clipped_end = 0
    match_splitted = re.split(r'(\d+)', sam[SAM_MID])
    match_splitted.pop(0) # re produces a weird '' as elem 0  ... 
    if match_splitted[1] == "S":
        soft_clipped_start = match_splitted[0]
    if match_splitted[-1] == "S":
        soft_clipped_end = match_splitted[-2]
    # take only the portion of the read that was aligned to the reference
    # (i.e. omit soft clippings at the start and at the end)
    last_pos = len(sam[SAM_READ]) - int(soft_clipped_end)
    r = sam[SAM_READ][int(soft_clipped_start):last_pos]
    r = list(r)

    # save the final reference and read sequences
    G = []
    R = []

    current_MID = sam[SAM_MID]
    exploded_MID = re.split(r'(\d+|[A-Za-z])', current_MID)
    for item in exploded_MID:
        if item == '':
            exploded_MID.remove(item)

    ### Do alignment read vs genome
    for i in range(int(len(exploded_MID)/2)):
        number = int(exploded_MID[i*2])
        type = exploded_MID[i*2 +1]

        # no need to process type "N" at this point
        # because currently the genomic seq anyways
        # only consists of "?"
        if type == "M":
            # read and reference were aligned
            for _ in itertools.repeat(None, number):
                G.append(g.pop(0))
                R.append(r.pop(0))                        
        elif type == "I":
            # nucleotides in read that are not there in the reference
            for _ in itertools.repeat(None, number):
                G.append("-")
                R.append(r.pop(0))
        elif type == "D":
            # nucleotides in reference that are not there in the read
            for _ in itertools.repeat(None, number):
                G.append(g.pop(0))
                R.append("-")


    mismatch = sam[SAM_MISMATCHES]
    mismatch = mismatch.split(":")[-1]
    mismatch = mismatch.replace("^", "")

    mismatch_splitted = re.split(r'(\d+|[A-Za-z])', mismatch)
    for item in mismatch_splitted:
        if item == '':
            mismatch_splitted.remove(item)

    pointer = 0

    for c in mismatch_splitted:

        if c == '':
            continue # Skip c and this iteration when this happens

        # Skip the next G[pointer] that is not a dash
        if c.isdigit():
            # perfect match here
            value = int(c)
            while value > 0:
                # skip gaps in the reference because insertions are
                # not present in the mismatch string
                if G[pointer] == '-':
                    pointer += 1
                else:
                    # go over a region of perfect matches
                    pointer += 1
                    value -= 1
        else:
            # again
            # skip gaps in the reference because insertions
            # are not present in the mismatch string
            while G[pointer] == '-':
                pointer += 1 
            # exchange question mark in the reference
            # two possible reasons: 
            # 1. mismatch; 2. deletion
            G.pop(pointer)
            G.insert(pointer, c)
            pointer+=1


    # for the MMID, track the positions of the mismatches and deletions
    # final MMID looks like:
    # 15MATIGDCDT12 (15 matches -> 1 Mismatch( A ref, T in read) -> 1 insertion (G in read) -> 1 deletion (C in reference) -> 1 deletion (T in reference) -> 12 matches
    match_count = 0
    MMID = ""


    g_index = 0
    r_index = 0

    for index in range(len(G)):
        cg = G[index]
        cr = R[index]
        if cg == "?":
            match_count += 1
        else:
            if match_count > 0:
                MMID += str(match_count)
                match_count = 0

            if cg == "-":
                MMID += "I"+cr

            elif cr == "-":
                MMID += "D"+cg

            else:
                MMID += "M"+cg+cr

        # note by RS 2016-06-01: this increment statement
        # should not effect the loop as far as I understood python
        # (was introduced not by me but by Manuel)
        index += 1

    if match_count > 0:
        MMID += str(match_count)
        match_count = 0

    index_space = sam[SAM_INDEX_SPACE]
    index_space = index_space.split(":")[-1]

    new_bed = list(bed)
    # exchange read weight by edit distance and
    # define the read weight later on
    new_bed.pop(4)
    new_bed.insert(4, index_space)
    new_bed.append(MMID)
    new_file_line = "\t".join(new_bed)
    return new_file_line

# Parses a string containing Matches (M) and InDels (I and D) and returns the number of positions represented by that string
def get_match_length( match_string):
    number_string = ""
    sequence_length = 0

    # the "r" in re.split(r...) means that the pattern is a raw string literal
    # so no escaping is done by python with the string
    # in the context of patterns this is helpful because it makes it 
    # obsolete to escape every backslash with another backslash
    match_splitted = re.split(r'(\d+)', match_string)
    match_splitted.pop(0) # re produces a weird '' as elem 0  ... 

    for i in range(int(len(match_splitted)/2)):
        value = match_splitted[i*2]
        type = match_splitted[i*2 + 1]

        # M -> match, D -> deletion in the read (hence, present in the reference)
        # N -> skipped in the read (hence, present in the reference; introns)
        if type == "M" or type == "D" or type == "N":
            sequence_length += int(value)
        else:
            pass

    return sequence_length

                
align()
