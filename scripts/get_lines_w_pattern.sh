#!/bin/bash

# Setting to strict mode
set -euo pipefail
IFS=$'\n\t'

#### FUNCTIONS

usage()
{
    echo "usage: get_lines_w_pattern.sh [[[-f file ] [-c column] [-p pattern] [-o outname]] | [-h]]"
}



#### MAIN
# Test whether all required inputs are present
if [[ $1 == -h ]] || [[ $# != 8 ]]
then
	usage
	exit
fi

# Get parameters
while [ $# -gt 0 ]; do
    case $1 in
        -f | --file )           shift
                                filename=$1
                                ;;
        -o | --out)             shift
                                outname=$1
                                ;;
        -c | --column )   	shift
				column=$1
                                ;;
        -p | --pattern )   	shift
				pattern=$1
                                ;;
        -h | --help )           usage
                                exit
                                ;;
        * )                     usage
                                exit 1
    esac
    shift
done

printf "\nAll lines with pattern \"%s\" in column %i are retrieved from %s and written to %s.\n" \
	${pattern}\
	${column}\
	${filename}\
	${outname}

# Get lines with matching pattern in requested column
awk -v col="${column}" -v pat="${pattern}" '$col == pat' ${filename} > ${outname}

printf "\nDONE!\n"
