#!/bin/sh -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y

# Time-stamp: <2015-05-16 20:50:29 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

## load perl module
module load perl-5.14


usage() {
    echo "usage: bedtobigwignormRPM.sh <LIB> <INDEXBASE> <ASSEMBLY>"
    echo "       Check that the 1bp window interval file matches"
    echo "       the mapped bed file."
    echo "where: <LIB> is a list of libraries to process"
    echo "       1bp library specific windows would be made."
    echo "This script does the job of making bedGraph files"
    echo "from normRPMdeciID bed files"
    echo "with the sum of the score for the read aligments"
    echo "and bigwig files."
}

# Minimal argument checking
if [ $# -lt 3 ]; then
    usage
    exit
fi


# Requirements
# bedops sort-bed, gnu awk, gnu sed, gnu coreutils, bedtools, UCSC Kent utilities

# Set variables
LIB=$1
INDEXBASE=$2
ASSEMBLY=$3

ASSEMBLYCHROMINFO=/data/OkamuraLab/local/processannot/${ASSEMBLY}.chromInfo.txt
NORMRPMID=mapcol${INDEXBASE}normRPMdeciID.bed
BWTYPE=normRPM

echo "Start"
date

echo "head check ${LIB}${NORMRPMID}"
head ${LIB}${NORMRPMID}
wc -l ${LIB}${NORMRPMID}


echo "make ${LIB} specific 1bp windows"
echo "do a bedtools merge on ${LIB}${NORMRPMID}"
bedtools merge -i ${LIB}${NORMRPMID} > ${LIB}${NORMRPMID}mergetemp

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}${NORMRPMID}mergetemp > ${LIB}${NORMRPMID}merge

rm ${LIB}${NORMRPMID}mergetemp

echo "head check ${LIB}${NORMRPMID}merge"
head ${LIB}${NORMRPMID}merge

echo "make window intervals of 1bp from the input bedfile"
bedtools makewindows -b ${LIB}${NORMRPMID}merge -w 1 > 1bp${INDEXBASE}.${LIB}specific.txt


echo "Do for ${INDEXBASE}"
echo "Check window intervals of 1bp, done" # do once
echo "head 1bp${INDEXBASE}.${LIB}specific.txt"
head 1bp${INDEXBASE}.${LIB}specific.txt

echo "Do for NormRPM"
echo "Split bed file with score to plus and minus files"
echo "plus strand"
awk '$6=="+"' ${LIB}${NORMRPMID} > ${LIB}${NORMRPMID}plustemp

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}${NORMRPMID}plustemp > ${LIB}${NORMRPMID}plus

rm ${LIB}${NORMRPMID}plustemp

echo "head check ${LIB}${NORMRPMID}plus"
head ${LIB}${NORMRPMID}plus

echo "minus strand"
awk '$6=="-"' ${LIB}${NORMRPMID} > ${LIB}${NORMRPMID}minustemp

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}${NORMRPMID}minustemp > ${LIB}${NORMRPMID}minus

rm ${LIB}${NORMRPMID}minustemp

echo "head check ${LIB}${NORMRPMID}minus"
head ${LIB}${NORMRPMID}minus

echo "Map to plus and minus files"
echo "add minus sign to value for minus strand bedGraph"

date

echo "Do for minus bedgraph"
bedtools map -o sum -c 5 -null 0 -a 1bp${INDEXBASE}.${LIB}specific.txt -b ${LIB}${NORMRPMID}minus | \
awk 'BEGIN {OFS=FS="\t"} {$4=sprintf("%.4f",$4)} 1' > ${LIB}${NORMRPMID}minus.bedGraph

echo "head check ${LIB}${NORMRPMID}minus.bedGraph"
head ${LIB}${NORMRPMID}minus.bedGraph
wc -l ${LIB}${NORMRPMID}minus.bedGraph

echo "Add minus sign to ${LIB}${NORMRPMID}minus.bedGraph"

awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,"-",$4}' ${LIB}${NORMRPMID}minus.bedGraph > ${LIB}${NORMRPMID}minus.bedGraphtemp01
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4$5}' ${LIB}${NORMRPMID}minus.bedGraphtemp01 > ${LIB}${NORMRPMID}minus.bedGraphtemp02

rm ${LIB}${NORMRPMID}minus.bedGraph
rm ${LIB}${NORMRPMID}minus.bedGraphtemp01

awk 'BEGIN {OFS=FS="\t"} {if ($4 == -0.0000) $4 = 0.0000; print}' ${LIB}${NORMRPMID}minus.bedGraphtemp02 > ${LIB}${NORMRPMID}minus.bedGraph.chr

rm ${LIB}${NORMRPMID}minus.bedGraphtemp02

echo "sort"

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}${NORMRPMID}minus.bedGraph.chr > ${LIB}${NORMRPMID}minus.bedGraph.chr.sort
mv ${LIB}${NORMRPMID}minus.bedGraph.chr.sort ${LIB}${NORMRPMID}minus.bedGraph

rm ${LIB}${NORMRPMID}minus.bedGraph.chr

date

echo "head check ${LIB}${NORMRPMID}minus.bedGraph"
head ${LIB}${NORMRPMID}minus.bedGraph

echo "check that the number of lines are correct"
wc -l ${LIB}${NORMRPMID}minus.bedGraph


echo "Do for plus bedgraph"
bedtools map -o sum -c 5 -null 0 -a 1bp${INDEXBASE}.${LIB}specific.txt -b ${LIB}${NORMRPMID}plus | \
awk 'BEGIN {OFS=FS="\t"} {$4=sprintf("%.4f",$4)} 1' > ${LIB}${NORMRPMID}plus.bedGraph

echo "sort"

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}${NORMRPMID}plus.bedGraph > ${LIB}${NORMRPMID}plus.bedGraph.sort
mv ${LIB}${NORMRPMID}plus.bedGraph.sort ${LIB}${NORMRPMID}plus.bedGraph

echo "head check ${LIB}${NORMRPMID}plus.bedGraph"
head ${LIB}${NORMRPMID}plus.bedGraph


echo "this is similar to bedgraph, with sum of read alignment score in bed file"
head ${LIB}${NORMRPMID}plus.bedGraph
head ${LIB}${NORMRPMID}minus.bedGraph


echo "convert bedgraph to bigwig files, check that bedgraph is sorted"
echo "use UCSC tool bedGraphToBigWig"

bedGraphToBigWig ${LIB}${NORMRPMID}plus.bedGraph ${ASSEMBLYCHROMINFO} ${LIB}.${ASSEMBLY}.${BWTYPE}plus.bw 
bedGraphToBigWig ${LIB}${NORMRPMID}minus.bedGraph ${ASSEMBLYCHROMINFO} ${LIB}.${ASSEMBLY}.${BWTYPE}minus.bw 

echo "remove unnecessary files"
rm ${LIB}${NORMRPMID}plus.bedGraph
rm ${LIB}${NORMRPMID}minus.bedGraph

rm ${LIB}${NORMRPMID}plus
rm ${LIB}${NORMRPMID}minus


rm 1bp${INDEXBASE}.${LIB}specific.txt
rm ${LIB}${NORMRPMID}merge

echo "Done"
date
