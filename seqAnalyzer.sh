#!/bin/bash

# usage seqAnalyzer.sh FastQreadfile DestinationFolder barcodeTextfile mappingFile 
# 
#
#
#

# BC : Barcode Split
# BR : Barcode Total Reads (high number)
# CU : Cutout
# CR : Cutout Total Reads
# CO : collapsed Sequences to fasta
# FR : Fasta Total Reads
# HF : BLAT Hits Filtered
# HR : Hits Total
# SC : Second collapse unique BLAT hits
# SR : Second collapse total Hits
# LN : Second collapse BLAT Reads
# LR : Second collapse total Reads
# MB : mask(loop) in barcode split files
##

#debugSwitches:
trimming=true
cutting=true
collapsing=true
library=true
blating=true
sorting=true
collapsing2=true
maskcheck=true
statistic=true


#default Value init
rThreshold=1000 # read threshold after barcode splitting
mmLoop=2
gapsLoop=2
keep="F"
itags="F"
loop="TAGTGAAGCCACAGATGTA"
mirSeq="GTATATTGCTGTTGACAGTGAGCG" #default 5'mirE
workdir=false
####################
while [[ $# > 0 ]];
do
    key="$1"
    case $key in
        -h|--help)
        ;;
        -k|--keepData)
        keep=true
        ;;
        -t|--threshold)
        rThreshold="$2"
        shift
        ;;
        -s|--missmatchLoop)
        mmLoop="$2"
        shift
        ;;
        -a|--gapsLoop)
        gapsLoop="$2"
        shift
        ;;
        -i|--itags)
        itags=true
        ;;
        -l|--loop)
        loop="$2"
        shift
        ;;
        -m|--mirCheck)
        mirSeq="$2"
        shift
        ;;
        -w|--workingDir)
        workdir="$2"
        shift
        ;;
        -j|--justSorting)
        keep=true
        trimming=true
        cutting=false
        collapsing=false
        library=false
        blating=false
        sorting=false
        collapsing2=false
        maskcheck=false
        statistic=false
        ;;
        *)
        break
        ;;
    esac
    shift
done

#### init variables
sourceFile="$1"
destination="$2"
barcodeFile=$3
mappingFile=$4
blatDestination="/home/mitschke/Data/Programms/BLAT/"
itaglength=7

adapterMask=""

echo -e "Source File  =\t${sourceFile}"
echo -e "Destination  =\t${destination}"
echo -e "Barcode File =\t${barcodeFile}"
echo -e "Library File =\t${mappingFile}"

absdest=$(readlink -e $(dirname "$destination"))
absdest="${absdest}/$(basename "${destination}")"
echo -e "Absolute d.  =\t${absdest}"
if [[ "$workdir" == false ]]; then
    #workdir="/tmp/seqAnalyzer_$$/"
    workdir="/Data/PipeTmp/seqAnalyzer_$$/"
    #workdir="/Data/PipeTmp/seqAnalyzer3_sumUp/"
    #workdir="/Data/PipeTmp/seqAnalyzer_25343/"
fi
mkdir -p $workdir
statFile=${workdir}Stats.txt
echo -e "#${1}\t${2}\t${3}\t${4}" > $statFile
#########
##Stats##

#tReads=$(expr $(cat "$sourceFile" | wc -l ) / 4) #| gawk '{if($2 !="total") print $2;}')
tReads=$(expr $(wc -l "$sourceFile" | gawk '{print $1}' ) / 4) #| gawk '{if($2 !="total") print $2;}')
#Total Reads
echo -e "TR\t$tReads" >> $statFile

scriptSource="/home/mitschke/Data/Scripts/shRNASeqAnalyzer/subScripts/"

## 1.###############
# SeqTrim
# options etc... -o outputfile -b barcodes -a adapter -s silentmode
if [[ "$trimming" == true ]]; then
    echo "Trimming ..."
    ${scriptSource}seqTrim.py -c "$barcodeFile" -d "${workdir}" -o "$1" -s "$sourceFile"
fi

## 2.###############
# HighBarcodes
# paste all Files with more than $rThreshold reads in a txt file, default 100000
# default now 0 
wc -l ${workdir}/*BC_????????.fastq | gawk -v t=${rThreshold} '$1 > t {print $0}' > "${workdir}_highBarcodes.txt"

## 3.###############
# SeqCutOut
# guide loop target cutout of highest barcodes
highBC=$(gawk ' {if($2 !="total") print $2;}' "${workdir}_highBarcodes.txt") 
if [[ "$cutting" == true ]]; then
    echo "Cutting ...${gapsLoop}!!"
    if [[ "$itags" == true ]]; then
        echo "Cutting for iTags!"
        ${scriptSource}seqCutOut.py -m "i" -g "$gapsLoop" -c "$barcodeFile" $highBC
    else
        echo "Cutting without iTags!"
        ${scriptSource}seqCutOut.py -g "$gapsLoop" -c "$barcodeFile" $highBC
    fi
fi

## 4.###############
# fastx_collapser
# collapse sequence reads
if [[ "$collapsing" == true ]]; then
    echo "Collapsing ..."
    for file in $(ls ${workdir}*_cutout.fasta);
    do
        fastx_collapser -i $file -o ${file%.*}_collapsed.fasta
    done


    ## 4a.###############
    # getItags.sh
    # set the iTags into the header
    if [[ "$itags" == true ]];
    then
        for cFile in $(ls ${workdir}*_collapsed.fasta);
        do
            echo -e "tagging .... ${cFile}"
            ${scriptSource}getiTags.sh -i "${cFile}" -o "${cFile%.fasta}_itags.fasta"
            echo -e "moving .... ${cFile}"
            mv "${cFile%.fasta}_itags.fasta" "$cFile"
        done
    fi
fi

## 4b.###############
# sortFasta
#split mappingFile into Pools
if [[ "$library" == true ]]; then
    echo "Library preparation ..."
    ${scriptSource}sortFasta.py -d $workdir "$mappingFile"
fi

## 5.###############
# BLAT 
if [[ "$blating" == true ]]; then
    col2blat=$(ls ${workdir}*_collapsed.fasta)

    for line in $col2blat; do
        barcode=${line%_cutout_collapsed.fasta}
        barcode=${barcode#*_BC_}
        poolfile=$(gawk -v b=$barcode '{if($2==b){ print $3; exit}}' "$barcodeFile")
        echo $poolfile
        if [[ $poolfile != "" ]]; then
            echo "$line Blat against $poolfile"
            ${blatDestination}blat -oneOff=1 -t=dna -q=dna -out=pslx -dots=5000 \
                 "${workdir}${poolfile}_pool.fasta" "$line" "${workdir}${barcode}_Pool(${poolfile})_BL.pslx"
        else
            echo "$line Blat agains All"
            ${blatDestination}blat -oneOff=1 -t=dna -q=dna -out=pslx -dots=5000 \
                 "${mappingFile}" "$line" "${workdir}${barcode}_Pool(All)_BL.pslx"
        fi
    done
fi

## 6.##############
#filter blat results
if [[ "$sorting" == true ]]; then
    echo "sorting ..."
    inputF=$(ls ${workdir}*_BL.pslx)
    ${scriptSource}sortBlatOutput.sh -n 10 -l "$mappingFile" $inputF
fi

## 7.#############
#collapse Blat results again and plots graphs
if [[ "$collapsing2" == true ]]; then
    echo "collapsing2 ..."
    inputP=$(ls ${workdir}*_HF.pslx)
    ${scriptSource}collapseBlat.py -s -c "$barcodeFile" -d "$workdir" $inputP 
fi

###############
###############
# Statistics! #


########## Stats File establishment

cat "${workdir}_highBarcodes.txt" | \
    gawk 'split($2,a,"BC_") split(a[2],b, ".") {if (b[1]!="") printf("BC\t%s\t%s\n", b[1], $1/4); else print "BR\t"$1/4;}'\
    >> $statFile

wc -l ${workdir}*_cutout.fasta | \
    gawk 'split($2,a,"BC_") split(a[2],b, "_cutout") {if (b[1]!="") print "CU\t"b[1]"\t"$1/2; else print "CR\t"$1/2;}'\
    >> $statFile

wc -l ${workdir}*_collapsed.fasta | \
    gawk 'split($2,a,"BC_") split(a[2],b, "_cutout") {if (b[1]!="") print "CO\t"b[1]"\t"$1/2; else print "FR\t"$1/2;}'\
    >> $statFile

##### BLAT Results 1. 
totalNum=0
for file in $(ls "${workdir}"*_HF.pslx);
do
    hits=$(wc -l $file | gawk '{print $1}')
    number=$(gawk 'BEGIN {n=0}; split($10,a,"-") {(n=n+a[2])}; END {print n}' $file) 
    barName=${file##*/}
    barName=${barName%_Pool*_HF.pslx}
    echo -e "HF\t${barName}\t${number}\t${hits}" >> $statFile
    let totalNum=$totalNum+$number
done
echo -e "HR\t$totalNum" >> $statFile

#wc -l ${workdir}*_HF.pslx | \
#    gawk 'split($2,a,"/") split(a[length(a)],b, "_Pool") {if (b[1]!="") print "HF\t"b[1]"\t"$1; else print "HR\t"$1;}'\
#    >> $statFile

##### BLAT Results 2.

wc -l ${workdir}*_SC.pslx | \
    gawk 'split($2,a,"/") split(a[length(a)],b, "_Pool") {if (b[1]!="total") print "SC\t"b[1]"\t"$1; else print "SR\t"$1;}'\
    >> $statFile

sndTotal=0
for sndFile in $(ls "${workdir}"*_SC.pslx); do
    readNumber=$(gawk 'BEGIN {n=0};split($10,a,"-") {n=n+a[2]}; END {print n}' "$sndFile")
    Code=${sndFile##*/}
    echo -e "LN\t${Code%_Pool*_SC.pslx}\t${readNumber}" >> $statFile
    #echo -e "LN\t${Code%_Pool*_SC.pslx}\t${readNumber}"

    let sndTotal=$sndTotal+$readNumber
done
echo -e "LR\t${sndTotal}" >> $statFile

#####
# check for shRNA mask (sense+loop+guide [22+19+22])
# and check for 5'mirE if set
# in outsorted 
# 1. in Barcode splitting
if [[ "$maskcheck" == true ]]; then
    echo -e "maskcheck"
    for outsorted in $(cat "${workdir}_highBarcodes.txt" | gawk '{if ($2!="total"){print $2}}'); do
        maskPos=$(gawk -v l=$loop '$0 ~ l' $outsorted | wc -l)
        mirNum=$(gawk -v m=$mirSeq '$0 ~ m' $outsorted | wc -l)
        code=${outsorted#*_BC_}
        code=${code%.fastq}
        #echo "${outsorted} for ${code} : ${maskPos}"
        echo -e "MB\t${code}\t${maskPos}"
        echo -e "MB\t${code}\t${maskPos}" >> $statFile
        echo -e "VB\t${code}\t${mirNum}"
        echo -e "VB\t${code}\t${mirNum}" >> $statFile
    done
fi

# 2. in Blat hits

##############
# Stat Maker #
##############
if [[ "$statistic" == true ]]; then
    if [[ "$itags" == true ]]; then
        ${scriptSource}statMaker.py -s "$statFile" -d "${workdir}" -e _pStat.b -c "$barcodeFile" -i
    else
        ${scriptSource}statMaker.py -s "$statFile" -d "${workdir}" -e _pStat.b -c "$barcodeFile"
    fi

fi

#${scriptSource}statMaker.py -s "$statFile" -d "${workdir}" -e _pStat.b -c "$barcodeFile"

#spike in controlling?
# number of iTags
# for file in $(ls *_Pool\(*\)_SC.pslx); do num=$(grep -e "hshRPA3.1281" "$file" | wc -l); echo -e "${num}"; done
# read number iTags Spike ins
# for file in $(ls *_Pool\(*\)_SC.pslx); do num=$(gawk 'BEGIN {n=0}; $0 ~ /hshRPA3.1281/ {split($10,a,"-") ;n=a[2]+n;}; END {print n}' $file); echo -e "${num}"; done


###############
# cleanup!
echo $absdest
mkdir -p $absdest
cp "${workdir}"*.png "${absdest}/"
cp "${workdir}"*_SC.pslx "${absdest}/"
cp "$statFile" "${absdest}/Stats.txt"
cp "${workdir}"*.tsv "${absdest}/"
cp "${workdir}"*.h5 "${absdest}/"

if [[ "$keep" != true ]];
then
    rm -r $workdir
fi


