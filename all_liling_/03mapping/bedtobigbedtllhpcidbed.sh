#!/bin/sh -login
#$ -cwd
#$ -V

## Join standard error and standard output into 1 file output
#$ -j y
# Time-stamp: <2014-01-18 00:00:11 liling>

usage() {
    echo "usage: bedtobigbedtllhpcidbed.sh <LIB> <INDEXBASE> <ASSEMBLY>"
    echo "where: <LIB> is the library name"
    echo "       <INDEXBASE> is the base name for the bowtie1 index"
    echo "       <ASSEMBLY> is the assembly name"
    echo "       written in each row provided in a tab-separated text file."
    echo "This script does the job of making bigbed files from <LIB>mapcol<INDEXBASE>ID.bed"
    echo "and removing the bed file to save space."
    echo "Location to chrominfo hardcoded in script."
}

# Minimal argument checking
if [ $# -lt 3 ]; then
    usage
    exit
fi

# Requirements
# bedops sort-bed, gnu awk, gnu sed, gnu coreutils, bedtools, UCSC Kent utilities

# Set variables
LIB=$1               #library name
INDEXBASE=$2         #chr1hsardna
ASSEMBLY=$3          #chr1hsardna

ASSEMBLYCHROMINFO=/data/OkamuraLab/local/processannot/${ASSEMBLY}.chromInfo.txt

INPUTBEDIDFILE=${LIB}mapcol${INDEXBASE}ID.bed

echo "Start with ${LIB}"
date

echo "print variables"
echo "<LIB>"
echo $1

echo "<INDEXBASE>"
echo $2

echo "<ASSEMBLY>"
echo $3


echo "head check ${INPUTBEDIDFILE}"
head ${INPUTBEDIDFILE}

echo "convert score column to 0 value"
awk 'BEGIN {OFS=FS="\t"} {print $1,$2,$3,$4,"0",$6}' ${INPUTBEDIDFILE} > ${INPUTBEDIDFILE}.col5tozero

head ${INPUTBEDIDFILE}.col5tozero

echo "Convert ${INPUTBEDIDFILE}.col5tozero to bigbed with UCSC Kent Utility bedToBigBed"
bedToBigBed ${INPUTBEDIDFILE}.col5tozero ${ASSEMBLYCHROMINFO} ${LIB}.${ASSEMBLY}.ID.bb

echo "remove unnecessary file"
echo "remove ${INPUTBEDIDFILE}.col5tozero"
rm ${INPUTBEDIDFILE}.col5tozero

echo "remove ${INPUTBEDIDFILE}"
rm ${INPUTBEDIDFILE}

echo "Done with ${LIB}"
date
