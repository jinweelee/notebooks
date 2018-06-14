#!/bin/bash -login
#$ -cwd
#$ -V
#$ -j y

# Time-stamp: <2015-12-10 16:11:33 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

#####################################################################################################################
# This script does the job of converting a bigbed file to bed, 
# normalize with features,
# then to convert to bigbed and bigwig files for visualization on UCSC GB.
# note that the identifier format is different with additional length information encoded.
# requirements:
# UCSC Kent Utilities, bedtools, bedops, gnu coreutils, gawk, gnu sed
#####################################################################################################################

usage() {
    echo "usage: normalizebigbedwithperthousandfeaturereadcounts.sh <LIB> <FEATURENORM> <FEATURENORMCOUNT> <COUNTTYPE> <ASSEMBLY>"
    echo "where: <LIB> is the library name"
    echo "       <FEATURENORM> is the feature to be normalized against"
    echo "       <FEATURENORMCOUNT> is the count value for the feature"
    echo "       <COUNTTYPE> is the type of count value file"
    echo "       <ASSEMBLY> is the genome assembly"
    echo "       written in each row provided by a tab-separated text file."
    echo "input: input file is bigbed file previously produced from mapping run."
    echo "This script does the job of making normalized bigbed and bigwig files."
    echo "Location to assembly chrominfo hardcoded in script."
    echo "Note that some parameter hardcoding done."
}

# Minimal argument checking
if [ $# -lt 5 ]; then
    usage
    exit
fi

## Set variables
LIB=$1
FEATURENORM=$2          #spikeinset01
FEATURENORMCOUNT=$3
COUNTTYPE=$4            #normRPKspikeinset01
ASSEMBLY=$5             #dm6
#INDEXBASE=dm6
REPORTPARAM=bowtiereportav0
LENGTHFILTER=l18t30nt

ASSEMBLYCHROMINFO=/data/OkamuraLab/local/processannot/${ASSEMBLY}.chromInfo.txt      ## TLL HPC
#ASSEMBLYCHROMINFO=/usr/local/processannot/${ASSEMBLY}.chromInfo.txt                 ## work computer

echo "Start"
date

echo "check locale setting"
locale

printf "\n"

echo "print variables"
echo "<LIB>"
echo $1
echo "<FEATURENORM>"
echo $2
echo "<FEATURENORMCOUNT>"
echo $3
echo "<COUNTTYPE>"
echo $4
echo "<ASSEMBLY>"
echo $5

printf "\n"

echo "This is running normalizebigbedwithfeature script from $PWD"

#####################################################################################################################
## bigbed to bed format

echo "use UCSC Kent Utility bigBedToBed to convert bigbed to bed"
bigBedToBed ${LIB}.${ASSEMBLY}.normReadID.bb ${LIB}.${ASSEMBLY}.normReadID.bed

echo "head check ${LIB}.${ASSEMBLY}.normReadID.bed"
echo "<CHR> <START> <STOP> <:${LIB}:sequence:countvalue:mapHit:lengthnt> <COUNTVALUE> <STRAND>"
head ${LIB}.${ASSEMBLY}.normReadID.bed

printf "\n"

echo "split identifier"
echo "<CHR> <START> <STOP> <lib> <sequence> <countvalue> <mapHit> <length> <STRAND>"
awk 'BEGIN {OFS=FS="\t"} {split($4,id,":"); print $1,$2,$3,id[2],id[3],id[4],id[5],length(id[3]),$6}' ${LIB}.${ASSEMBLY}.normReadID.bed \
> ${LIB}.${ASSEMBLY}.normReadID.txt

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.normReadID.txt"
echo "<CHR> <START> <STOP> <lib> <sequence> <countvalue> <mapHit> <length> <STRAND>"
head ${LIB}.${ASSEMBLY}.normReadID.txt

printf "\n"

#echo "remove ${LIB}.${ASSEMBLY}.normReadID.bed"
#rm ${LIB}.${ASSEMBLY}.normReadID.bed

printf "\n"

#####################################################################################################################
## normalize reads per thousand ${FEATURENORM} read counts

echo "normalize reads per thousand ${FEATURENORM} read counts"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,$5,sprintf("%.4f",(($6/'$FEATURENORMCOUNT')*1000)),$7,$8,$9}' ${LIB}.${ASSEMBLY}.normReadID.txt \
> ${LIB}.${ASSEMBLY}.${COUNTTYPE}.txt

echo "head check ${LIB}.${ASSEMBLY}.${COUNTTYPE}.txt"
head ${LIB}.${ASSEMBLY}.${COUNTTYPE}.txt

printf "\n"

echo "format to bed file"
echo "format identifier"
echo "<CHR> <START> <STOP> <:${LIB}:sequence:countvalue:mapHit:lengthnt> <COUNTVALUE> <STRAND>"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,":"$4":"$5":"$6":"$7":"$8"nt",$6,$9}' ${LIB}.${ASSEMBLY}.${COUNTTYPE}.txt > ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed

printf "\n"

#####################################################################################################################
## bed to bigbed

echo "make bigbed files"
echo "convert score column to 0 value"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,"0",$6}' ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed.col5tozero

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed.col5tozero"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed.col5tozero

printf "\n"

echo "convert ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed.col5tozero to bigbed with UCSC Kent Utility bedToBigBed"
bedToBigBed ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed.col5tozero ${ASSEMBLYCHROMINFO} ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bb

printf "\n"

echo "remove unnecessary file"
echo "remove ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed.col5tozero"
rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed.col5tozero

printf "\n"

#####################################################################################################################
## bed to bigwig

echo "make bigwig files"
echo "make ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed specific 1bp windows"
echo "do a bedtools merge on ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed"
bedtools merge -i ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmergetemp

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmergetemp \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmerge

rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmergetemp

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmerge"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmerge

printf "\n"

echo "make window intervals of 1bp from the input bedfile"
bedtools makewindows -b ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmerge -w 1 \
> 1bp.${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed

echo "head 1bp.${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed"
head 1bp.${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed

echo "remove unnecessary files"
rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedmerge

printf "\n"

echo "Split bed file with score to plus and minus files"
echo "plus strand"
awk '$6=="+"' ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplustemp

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplustemp \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus

rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplustemp

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus
wc -l ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus

printf "\n"

echo "minus strand"
awk '$6=="-"' ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminustemp

sort-bed --max-mem 4G --tmpdir $PWD ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminustemp \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus

rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminustemp

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus
wc -l ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus

printf "\n"
#echo "remove unnecessary files"
#rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed

printf "\n"

echo "Map to plus and minus files"
echo "do for plus bedgraph"
bedtools map -o sum -c 5 -null 0 -a 1bp.${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed -b ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus | \
    awk 'BEGIN {OFS=FS="\t"} {$4=sprintf("%.4f",$4)} 1' > ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph

echo "sort bed"
sort-bed --max-mem 4G --tmpdir $PWD ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph > ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph.sort
mv ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph.sort ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph

echo "remove unnecessary files"
rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus

printf "\n"

echo "do for minus bedgraph"
echo "add minus sign to value for minus strand bedGraph"
bedtools map -o sum -c 5 -null 0 -a 1bp.${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed -b ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus | \
    awk 'BEGIN {OFS=FS="\t"} {$4=sprintf("%.4f",$4)} 1' > ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph
wc -l ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph

echo "remove unnecessary files"
rm 1bp.${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bed

printf "\n"

echo "add minus sign to ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,"-",$4}' ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraphtemp01

awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4$5}' ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraphtemp01 \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraphtemp02

rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph
rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraphtemp01

awk 'BEGIN {OFS=FS="\t"} {if ($4 == -0.0000) $4 = 0.0000; print}' ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraphtemp02 \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph.chr

rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraphtemp02

echo "sort bed"
sort-bed --max-mem 4G --tmpdir $PWD ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph.chr \
> ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph.chr.sort

mv ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph.chr.sort ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph

rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph.chr

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph

printf "\n"

echo "check that the number of lines are correct"
wc -l ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph

printf "\n"

echo "this is similar to bedgraph, with sum of read alignment score in bed file"
echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph

printf "\n"

echo "head check ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph"
head ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph

printf "\n"

echo "remove unnecessary files"
rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus

printf "\n"

#####################################################################################################################
## bedgraph to bigwig format

echo "convert bedgraph to bigwig files, check that bedgraph is sorted"
echo "use UCSC tool bedGraphToBigWig"
bedGraphToBigWig ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph ${ASSEMBLYCHROMINFO} ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}plus.bw 
bedGraphToBigWig ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph ${ASSEMBLYCHROMINFO} ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}minus.bw     

printf "\n"

echo "remove unnecessary files"
rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedplus.bedGraph
rm ${LIB}.${ASSEMBLY}.${REPORTPARAM}.${LENGTHFILTER}.${COUNTTYPE}.bedminus.bedGraph

#####################################################################################################################


echo "Done"
date
