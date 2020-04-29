#!/usr/bin/env python
import sys
import re
import gzip
from argparse import ArgumentParser, RawTextHelpFormatter


### Created: Mar 5, 2019
### Author: Paula Iborra 
### Company: Zavolan Group, Biozentrum, University of Basel


### ARGUMENTS ###

parser = ArgumentParser(
    description="Script to filter FASTA files"
    )
parser.add_argument(
    '-v','--version', 
    action='version', 
    version='%(prog)s 1.0',
    help="Show program's version number and exit"
    )
parser.add_argument(
    '--trim', 
    help="Character's used to trim the ID. Remove anything that follows the character's. Write \\ infront of \'.\' and \'-\'(i.e trim=\"$\\.\\-|_\").  Default: first white space",
    type=str, 
    nargs='?', 
    default=""
    )
parser.add_argument(
    '--idlist', 
    help="Generate text file with the sequences IDs. One ID per line."
    )
parser.add_argument(
    '-f','--filter',
    help="Input ID list. Filter IDs and sequences from FASTA file with the mode selected. Filter file must contain ONE ID per line",
    )
parser.add_argument(
    '-m', '--mode',
    help="Type of filtering fasta file: keep (k) or discard (d) IDs contained in the ID list file.",
    choices=('k', 'd')
    )
parser.add_argument(
    '-r','--remove',
    help="Remove sequences from FASTA file longer than specified length.",
    type=int
    )
parser.add_argument(
    '-i','--input', 
    required=True, 
    help="Input FASTA file", 
    type=str
    )
parser.add_argument(
    '-o','--output', 
    help="Output FASTA file"
    )
    
args = parser.parse_args()

if args.filter and not args.mode:
    sys.exit("ERROR! Mode argument required when using filter option. (--mode, -m). See --help option.")

### PARSE FASTA FILE ###

class Seq:
    def __init__(self):
        self.id=""
    def __init__(self):
        self.seq=""
    def __init__(self):
        self.features=""

# open files 
if args.input.endswith('.gz'):
    f = gzip.open(args.input, 'rt')
else:
    f = open(args.input)

record=[] #list of records 
nrec=-1
inseq=0

# parse fasta file
sys.stdout.write("Parsing FASTA file...")
for line in f:
    if re.match(r'^>' , line):
        nrec+=1
        record.append(Seq())

        # define id of the record
        if not args.trim:
            mobj=re.match(r'^>(\S*)(.*)', line) 
        else:
            mobj=re.match(r'^>([^%s]*)(.*)'%args.trim , line)

        # add id and features
        if (mobj):
            record[nrec].id=mobj.group(1)
            record[nrec].features=mobj.group(2)
        inseq=0
    else :
        if inseq==0 :
            inseq=1
            record[nrec].seq = line   
        else:
            cstring=record[nrec].seq+line
            record[nrec].seq = cstring
sys.stdout.write("DONE\n")

## ID FILTER LIST ##

if (args.filter):
    sys.stdout.write("Filtering FASTA file...")
    id_filter=[line.rstrip('\n') for line in open(args.filter)]
    sys.stdout.write("DONE\n")


## OUTPUT FASTA FILE ##

if (args.output):
    sys.stdout.write("Writing FASTA file...")
    with open(args.output, 'w') as output:
        if (args.filter) and args.mode == 'k':
            if (args.remove):
                for x in range(0,nrec+1):
                    if record[x].id in id_filter and (len(record[x].seq)-1 <= args.remove):
                        output.write(">%s\n%s"%(record[x].id, record[x].seq))
            else:
                for x in range(0,nrec+1):
                    if record[x].id in id_filter:
                        output.write(">%s\n%s"%(record[x].id, record[x].seq))
        elif (args.filter) and args.mode == 'd':
            if (args.remove):
                for x in range(0,nrec+1):
                    if record[x].id not in id_filter and (len(record[x].seq)-1 <= args.remove):
                        output.write(">%s\n%s"%(record[x].id, record[x].seq))
            else:
                for x in range(0,nrec+1):
                    if record[x].id not in id_filter: 
                        output.write(">%s\n%s"%(record[x].id, record[x].seq))
        else:
            if (args.remove):
                for x in range(0,nrec+1):
                    if (len(record[x].seq)-1 <= args.remove):
                        output.write(">%s\n%s"%(record[x].id, record[x].seq))
            else:
                for x in range(0,nrec+1):
                    output.write(">%s\n%s"%(record[x].id, record[x].seq))
    output.close()
    sys.stdout.write("DONE\n")


## OUTPUT LIST IDs ##    

idlist=[]
if (args.idlist):
    sys.stdout.write("Creating IDs list from FASTA file...")
    fasta = open(args.output, 'r')
    with open(args.idlist, 'w') as id_list:
        for line in fasta:
            if line.startswith('>'):
                idlist.append(line[1:])
        idlist.sort()
        id_list.write(''.join(idlist))
        id_list.close()
        sys.stdout.write("DONE\n")
