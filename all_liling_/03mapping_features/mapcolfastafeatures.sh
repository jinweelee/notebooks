#!/bin/bash -login
#$ -cwd
#$ -V
#$ -pe smp 2
#$ -j y
# Time-stamp: <2015-12-03 16:31:25 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

## load perl module
module load perl-5.14

usage () {
    echo "usage: mapcolfastafeatures.sh <LIB> <MINLENGTH> <MAXLENGTH> <MAPINDEX> <INDEXBASE> <NTMISMATCH> <ALIGNREPORT> <FXARTIFACTFILTER>"
    echo "where: <LIB> is the library name"
    echo "       <MINLENGTH> is the minimum read length"
    echo "       <MAXLENGTH> is the maximum read length"
    echo "       <MAPINDEX> is the bowtie1 index to map"
    echo "       <INDEXBASE> is the base name for the bowtie1 index"
    echo "       <NTMISMATCH> is the number of nucleotide mismatches allowed"
    echo "       <ALIGNREPORT> is the number of alignments to report"
    echo "       <FXARTIFACTFILTER> is the fastx artifact filtered file to remove reads with all but three identical bases"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of mapping collapsed fasta files"
    echo "to bowtie1 index with a list of variables."
    echo "Location to fasta file and bowtie indexes hardcoded in script."
}

# Minimal argument checking
if [ $# -lt 7 ]; then
    usage
    exit
fi

# Set variables
LIB=$1
MINLENGTH=$2  #18
MAXLENGTH=$3  #30
MAPINDEX=$4   #mm10onlychr
INDEXBASE=$5  #mm10onlychr
NTMISMATCH=$6 #0
ALIGNREPORT=$7 #a
FXARTIFACTFILTER=$8 #fxaf

LENGTHFILTERFASTA=l${MINLENGTH}t${MAXLENGTH}nt

BOWTIE_INDEXES=/data/OkamuraLab/local/bowtieindexes/

INPUTMAPFASTAFILE=/data/OkamuraLab/processedlibraries/clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colidfasta/${LIB}clip${LENGTHFILTERFASTA}${FXARTIFACTFILTER}colid.fasta

#NSLOTS=1   # see pe smp above

echo "Start with ${LIB}"
date

echo "print variables"
echo "<LIB>"
echo $1
echo "<MINLENGTH>"
echo $2
echo "<MAXLENGTH>"
echo $3
echo "<MAPINDEX>"
echo $4
echo "<INDEXBASE>"
echo $5
echo "<NTMISMATCH>"
echo $6
echo "<ALIGNREPORT>"
echo $7
echo "<FXARTIFACTFILTER>"
echo $8

echo "check locale setting"
locale

echo "This is running mapcolfastafeatures script from $PWD"

printf "\n"

echo "head check ${INPUTMAPFASTAFILE} file"
head ${INPUTMAPFASTAFILE}

printf "\n"

echo "running bowtie, ${MAPINDEX} index, fasta input, \
report ${NTMISMATCH} mismatches, ${ALIGNREPORT} alignments, bowtie output format"
bowtie -f --un ${LIB}unmapcol${INDEXBASE}.fasta -${ALIGNREPORT} -v ${NTMISMATCH} -p ${NSLOTS} -t ${MAPINDEX} ${INPUTMAPFASTAFILE} \
> ${LIB}mapcol${INDEXBASE}.bowtie.txt

echo "parameters used"
echo "bowtie -f --un ${LIB}unmapcol${INDEXBASE}.fasta -${ALIGNREPORT} -v ${NTMISMATCH} -p ${NSLOTS} -t ${MAPINDEX}"

echo "done with mapping"

echo "head check ${LIB}mapcol${INDEXBASE}.bowtie.txt"
head ${LIB}mapcol${INDEXBASE}.bowtie.txt

echo "count how many reads were mapped for ${LIB}mapcol${INDEXBASE}.bowtie.txt"
MAP=$(awk '{OFS=FS="\t"} {a[$1]++} END {for (b in a) {print b}}' ${LIB}mapcol${INDEXBASE}.bowtie.txt | \
awk 'BEGIN {OFS=FS="\t"} {split($1,idm,"_"); print idm[1],idm[2],idm[3]}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $MAP

echo "count how many unmapped reads in ${LIB}unmapcol${INDEXBASE}.fasta"
UNMAP=$(fasta_formatter -i ${LIB}unmapcol${INDEXBASE}.fasta -t | \
awk 'BEGIN {OFS=FS="\t"} {split($1,idu,"_"); print idu[1],idu[2],idu[3]}' | awk '{SUM+=$3} END {if (SUM > 0) print SUM; else print "0"}')

echo $UNMAP

echo "sum of map and unmap reads"
SUMREAD=$(($MAP+$UNMAP))
echo $SUMREAD

echo "percentage of mapped reads with ${NTMISMATCH} mismatches"
PER=$(awk 'BEGIN {printf ("%0.2f", '$MAP' / '$SUMREAD' * 100)}')
echo $PER

printf "\n"

echo "This done running mapcolfastafeatures script from $PWD"

echo "Done with ${LIB}"
date
