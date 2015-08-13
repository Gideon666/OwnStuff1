#!/bin/bash

# get the Itags from the end of the sequence and put them into the header
# it just takes the last seven chars
# itag is deleted from the sequence
# usage:
#
# init
seq="CTAGAATT"

while [[ $# > 0 ]];
do
    key="$1"
    case $key in
        -h|--help)
        ;;
        -o|--output)
        outputFile="$2"
        shift
        ;;
        -i|--input)
        iFile="$2"
        shift
        ;;
        *)
        break
        ;;
    esac
    shift
done
echo "Input ${iFile}, Output $outputFile"
#for file in "$@"; do
    #outputFile="${file%_cutout_collapsed.fasta}"
    #outputFile="${outputFile}_itags_cutout_collapsed.fasta"
    gawk '$0 ~ /^>/ {a=$0; next}; {print a"-"substr($0, length($0)-6)"\n"substr($0,0,length($0)-7)}' "$iFile" >> "$outputFile"
#done
