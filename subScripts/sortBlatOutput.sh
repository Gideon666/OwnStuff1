#!/bin/bash
#gawk '($9 == "+") && (($1 + $4) == $11) {print $0}' ACGTACGT.output > ACGTACGT_hits.output
# usage sortBlatOutput.sh [options] -l databaseFile FileList


docline="##########################################################\n
usage: [options] blatresult(pslx file)\n
    \t-h, --help:\n
    \t\tthis help page.\n
    \t-l, --library:
    \t\t"




allowedN=10 # in loop
allowedMM=0
allowedGaps=0
allowedGBases=0
allowedHitLen=63
keep="False"
#### option read in
while [[ $# > 0 ]]
do
    key="$1"
        
    case $key in
        -h|--help)
        echo -e "$docline"
        exit 
        ;;
        -l|--library)
        DatabaseFile="$2"
        DatabaseEntrys=$(wc -l "$DatabaseFile")
        DatabaseEntrys=${DatabaseEntrys%% *}
        shift
        ;;
        -n|--allowedN)
        allowedN=$2
        shift
        ;;
        -m|--allowedM)
        allowedMM=$2
        shift
        ;;
        -a|--allowedHitLen)
        allowedHitLen=$2
        shift
        ;;
        -k|--keep)
        keep="True"
        ;;
        *)
        break
        ;;
    esac
    shift
done
#### noitpo

#entrys=expr $DatabaseEntrys / 2
#echo $DatabaseEntrys
#OLDIFS=$IFS
#IFS='.'
counter=0
barcodes=(AAGGTTCC GATCGATC CCAAGGTT TTCCAAGG GGTTCCAA GGAAAACC CCTTTTGG GGCCTTTT CCGGAAAA TTTTGGCC AAAACCGG ATCGATCG GTCAGTCA ACGTACGT GTACGTAC ATCGGCTA TAGCTTGT CGATGTTT GCCAATGT ACAGTGGT ATCACGTT GATCAGCG CAGATCTG TTAGGCAT GGCTACAG GGCTACAG CTTGTACT ACTTGATG TGACCACT TGGTTGTT TTCGCACC GCGCATAT TCTAGCTA AACCTGTG CACAGTGA GTGACACT TCGCTCAG TCAGATCT AATTGGCC GGGTTCCA CCCAAGGA)

# PSLX header
# 0      1      2       3       4       5       6       7       8       9               10      11      12      13              14      15      16      17      18              19      20  
# match  mis-   rep.    N's     Q gap   Q gap   T gap   T gap   strand  Q               Q       Q       Q       T               T       T       T       block   blockSizes      qStarts tStarts
#        match  match           count   bases   count   bases           name            size    start   end     name            size    start   end     count

echo $@

declare -A collapse
for bName in $@
do  
    #if [ $bName == 'output' ]
    #then
    #    continue
    #fi
    #if [ $bName in $barcodes ]
    #then
    #    poolNumber=$barcodes.index($bName)
    #else
    #then
    #    continue
    #fi
    #IFS=$OLDIFS
    #echo $bName
    fileholder[counter]="${bName}_hits.output"
    suffixHit="_HF.pslx"
    suffixNoHit="_BLATOUT.pslx"
    #gawk -v n=$allowedN -v m=$allowedMM \
    #    'split($14,a,":") split($22,b,"n") {if(($9 == "+") && (($1 + $4) == 63) && (($4 <= (n-1)) || (($4 == n) && (23<=length(b[1])) && length(b[1])<=42)) && ($2 <= m)) print $0;}' "${bName}" > "${bName%_BL.*}_HF.pslx"
    hitfile="${bName%_BL.*}_HF.pslx"
    outfile="${bName%_BL.*}_BLATOUT.pslx"
    gawk -v n=$allowedN -v m=$allowedMM -v ifile="$hitfile" -v ofile="$outfile" -v g="$allowedGaps" -v hl="$allowedHitLen"\
        '{if(($9 == "+") && (($1 + $4) >= hl) && (($4 <= n) && ($2 <= m) && ($7 <= g))) {print $0 > ifile} else {print $0 > ofile}}'\
       "${bName}"
    #gawk -v n=$allowedN -v m=$allowedMM '{if(not(($9 == "+") && (($1 + $4) == 63) && (($4 <= n) && ($2 <= m)))) print $0;}' "${bName}" > "${bName%_BL.*}_BLATOUT.pslx"
    perfecthits=$(cat ${bName%_BL.*}_HF.pslx | gawk '{print $14}' | sort -u | wc -l)
    #for line in $(sort -k 14 ${bName}_hits.output)
    #do
    #    linecounter=0
    #    item=$(gawk 'split($10,a,"-") {print $14, a[2]}' $line)
    #done
    ((counter++))
    echo ${bName%.*} $perfecthits
#   IFS='.'
done
#IFS=$OLDIFS
echo "finished sorting ${bName}..."
#echo ${fileholder[@]}
