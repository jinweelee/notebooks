#!/bin/sh -login
#$ -cwd
#$ -V
#$ -pe smp 2

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2015-05-16 20:48:43 liling>

## exit when there is an error
set -e

## set locale as C
export LC_ALL="C"

## load perl module
module load perl-5.14

usage() {
    echo "usage: getlistformapcolfastabwbbnormRPM.sh <DATAPROCESSLIST>"
    echo "where: <DATAPROCESSLIST> is a text file with each row for a library"
    echo "       to process with <LIB> <MINLENGTH> <MAXLENGTH> <MAPINDEX> <INDEXBASE> <NTMISMATCH> <ALIGNREPORT> <BWBBTYPE> <ASSEMBLY> <FXARTIFACTFILTER>"
    echo "       in which"
    echo "       <LIB> is the library name"
    echo "       <MINLENGTH> is the minimum read length"
    echo "       <MAXLENGTH> is the maximum read length"
    echo "       <MAPINDEX> is the bowtie1 index to map"
    echo "       <INDEXBASE> is the base name for the bowtie1 index"
    echo "       <NTMISMATCH> is the number of nucleotide mismatches allowed"
    echo "       <ALIGNREPORT> is the number of alignments to report"
    echo "       <BWBBTYPE> is the type of bigwig and bigbed file to make"
    echo "       <ASSEMBLY> is the assembly name"
    echo "       <FXARTIFACTFILTER> is the fastx artifact filtered file to remove reads with all but three identical bases"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of fetching the list of variables to pass to"
    echo "the script mapcolfasta.sh,"
    echo "and to make normRPM bigwig and bigbed files."
}


# Minimal argument checking

if [ $# -lt 1 ]; then
    usage
    exit
fi

# Set variable for input file
DATAPROCESSLIST=$1

echo "Start"
date

while read -r LINE
do
    LIB=$(echo $LINE | awk '{print $1}')
    MINLENGTH=$(echo $LINE | awk '{print $2}')
    MAXLENGTH=$(echo $LINE | awk '{print $3}')
    MAPINDEX=$(echo $LINE | awk '{print $4}')
    INDEXBASE=$(echo $LINE | awk '{print $5}')
    NTMISMATCH=$(echo $LINE | awk '{print $6}')
    ALIGNREPORT=$(echo $LINE | awk '{print $7}')
    BWBBTYPE=$(echo $LINE | awk '{print $8}')
    ASSEMBLY=$(echo $LINE | awk '{print $9}')
    FXARTIFACTFILTER=$(echo $LINE | awk '{print $10}')
    echo "variables to process"
    echo $LIB $MINLENGTH $MAXLENGTH $MAPINDEX $INDEXBASE $NTMISMATCH $ALIGNREPORT $BWBBTYPE $ASSEMBLY $FXARTIFACTFILTER
    echo "run mapping bowtie1"
    mapcolfasta.sh $LIB $MINLENGTH $MAXLENGTH $MAPINDEX $INDEXBASE $NTMISMATCH $ALIGNREPORT $FXARTIFACTFILTER
    echo "run bed to bedgraph to bigwig for normRPM"
    bedtobigwignormRPM.sh $LIB $INDEXBASE $ASSEMBLY
    echo "run bed to bigbed for <LIB>mapcol<INDEXBASE>ID.bed and remove bed file"
    bedtobigbedtllhpcidbed.sh $LIB $INDEXBASE $ASSEMBLY
    echo "run bed to bigbed for <LIB>mapcol<INDEXBASE>normRPMdeciID.bed and remove bed file"
    bedtobigbedtllhpcnormRPMidbed.sh $LIB $INDEXBASE $ASSEMBLY
    echo "run bed to bigbed for <LIB>mapcol<INDEXBASE>normReadID.bed and remove bed file"
    bedtobigbedtllhpcnormReadidbed.sh $LIB $INDEXBASE $ASSEMBLY
done < $1

echo "Done"
date
